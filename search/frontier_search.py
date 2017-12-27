#! /usr/bin/env python
import argparse
import os
import sys
from copy import copy

import sourmash_lib
from sourmash_lib import MinHash, signature
from sourmash_lib.sourmash_args import load_query_signature
from typing import Dict, List, Set, Union, Tuple

from .memoize import memoize
from . import search_utils
from .search_utils import get_minhashdb_name, load_dag, load_minhash
from spacegraphcats.logging import log


class NoContainment(ValueError):
    pass


def find_shadow(nodes: List[int], dag: Dict[int, List[int]]) -> Set[int]:
    shadow = set()  # type: Set[int]

    seen_nodes = set()  # type: Set[int]

    def add_to_shadow(node_id: int):
        if node_id in seen_nodes:
            return
        seen_nodes.add(node_id)

        children_ids = dag[node_id]

        if len(children_ids) == 0:
            shadow.add(node_id)
        else:
            for child in children_ids:
                add_to_shadow(child)

    for node in nodes:
        add_to_shadow(node)

    return shadow

def compute_overhead(node_minhash: MinHash, query_minhash: MinHash) -> float:
    """ Compute the relative overhead of minhashes.
    That is, the number of minhashes that are also in the query divided by the number of minhashes. """
    node_length = len(node_minhash.get_mins())
    return (node_length - node_minhash.count_common(query_minhash)) / node_length


def frontier_search(query_sig, top_node_id: int, dag, minhash_db: Union[str, search_utils.MinhashSqlDB], max_overhead: float, bf, vardb, include_empty = False, use_purgatory = False):
    """
        include_empty: Whether to include nodes with no minhases in the frontier
        use_purgatory: put leaf nodes with large overhead into a purgatory and consider them after we have processes the whole frontier
    """
    # expand the frontier where the child nodes together have more than x% overhead

    if not query_sig.minhash.scaled:
        raise Exception("query signature must be created with --scaled")

    # nodes and minhases that are in the frontier
    frontier = []  # type: List[int]
    frontier_minhash = None

    # nodes and minhashes for leaves with large overhead
    purgatory = []  # type: List[Tuple[float, int, MinHash]]

    seen_nodes = set()  # type: Set[int]

    num_leaves = 0
    num_empty = 0

    @memoize
    def var_in_bf(node_id):
        var_mh = load_minhash(node_id, vardb)
        if var_mh:
            mins = set(var_mh.get_mins())
            for hashval in mins:
                if bf.get(hashval):
                    return True
        return False

    @memoize
    def var_in_bf_decide(node_id):
        var_mh = load_minhash(node_id, vardb)
        if not var_mh:                    # nothing to decide upon, nokeep
            return False

        MAX_MINS=10                       # @CTB don't hardcode, maaan
        mins = set(var_mh.get_mins())
        sum_in = 0
        for hashval in mins:
            if bf.get(hashval):
                sum_in += 1

        frac = sum_in / len(mins)
        oh = 1 - frac

        is_full = False
        if len(mins) == MAX_MINS:
            is_full = True

        return is_full, frac, oh

    @memoize
    def load_and_downsample_minhash(node_id: int) -> MinHash:
        minhash = load_minhash(node_id, minhash_db)
        if minhash is None:
            return None
        max_scaled = max(query_sig.minhash.scaled, minhash.scaled)
        return minhash.downsample_scaled(max_scaled)

    @memoize
    def get_query_minhash(scaled: int) -> MinHash:
        return query_sig.minhash.downsample_scaled(scaled)

    @memoize
    def node_overhead(minhash: MinHash) -> float:
        query_mh = get_query_minhash(minhash.scaled)
        return compute_overhead(minhash, query_mh)

    @memoize
    def node_containment(minhash: MinHash) -> float:
        query_mh = get_query_minhash(minhash.scaled)
        return query_mh.contained_by(minhash)

    def add_node(node_id: int, minhash: MinHash):
        nonlocal frontier_minhash
#        if node_id == top_node_id:      # this should never happen
#            raise Exception             # ...unless the match is whole graph!

        frontier.append(node_id)
        if minhash:
            if frontier_minhash:
                minhash = minhash.downsample_scaled(frontier_minhash.scaled)
                frontier_minhash.merge(minhash)
            else:
                frontier_minhash = copy(minhash)

    def add_to_frontier(node_id: int):
        """
        Add a node or its children to the frontier.
        """
        if node_id in seen_nodes:
            # we have already added this node or its children to the frontier
            return
        else:
            seen_nodes.add(node_id)

        minhash = load_and_downsample_minhash(node_id)

        children_ids = dag[node_id]
        if len(children_ids) == 0:
            # the node is a leaf so let's add it to the frontier
            nonlocal num_leaves
            if use_purgatory:
                # nodes with high overhead are added to the purgatory for later evaluation whether
                # they are actually necessary to cover the query
                overhead = node_overhead(minhash)
                if overhead > max_overhead:
                    purgatory.append((overhead, node_id, minhash))
                else:
                    add_node(node_id, minhash)
                    num_leaves += 1
            else:
                add_node(node_id, minhash)
                num_leaves += 1
            return

        # check whether the node has more than x% overhead

        containment = node_containment(minhash)

        if containment == 0:
            # ignore nodes that have no containment
            # sanity check, should not get here
            raise NoContainment("No containment")

        overhead = node_overhead(minhash)

        # print("{} Overhead {}".format(node_id, overhead))

        if overhead > max_overhead:
            # greedily find the minimum number of children such that they together contain everything the node contains

            overheads = []

            for child_id in children_ids:
                child_mh = load_and_downsample_minhash(child_id)
                if not child_mh:
                    if include_empty:
                        nonlocal num_empty
                        # always add nodes without minhashes to frontier
                        add_node(child_id, None)
                        num_empty += 1
                    continue

                # ignore children without any containment
                if node_containment(child_mh) == 0:
                    continue

                child_oh = node_overhead(child_mh)
                overheads.append((child_oh, child_id, child_mh))

            # look at nodes one by one, starting with the smallest overhead and try to cover what the parent covered
            overheads.sort()

            _, first_node, first_mh = overheads[0]

            add_to_frontier(first_node)

            # Get a list of the minhases that need to be covered because they are in the query and in the parent node.
            query_mh = get_query_minhash(first_mh.scaled)
            required_query_minhashes = set(minhash.get_mins()).intersection(query_mh.get_mins()).difference(first_mh.get_mins())

            # Repeatedly remove hashes that are covered by nodes and add nodes to frontier until no more minhashes need to be covered.    
            for _, child_id, child_mh in overheads[1:]:
                if len(required_query_minhashes) == 0:
                    # early termination, all children already cover the node so we can stop
                    return

                before = len(required_query_minhashes)
                required_query_minhashes -= set(child_mh.get_mins())
                after = len(required_query_minhashes)
                if before > after:
                    add_to_frontier(child_id)

            if len(required_query_minhashes) > 0:
                raise Exception('Children cannot cover node.')

        else:
            # low overhead node gets added to the frontier
            add_node(node_id, minhash)

    def child_is_interesting(node_id):
        is_int = False

        children_ids = dag[node_id]
        if not children_ids:
            if var_in_bf(node_id):
                return True
            return False

        for child_id in children_ids:
            if child_is_interesting(child_id):
                return True

    def add_to_frontier2(node_id):
        if node_id in seen_nodes:
            return
            
        seen_nodes.add(node_id)

        minhash = load_and_downsample_minhash(node_id)

        if minhash:                     # non-empty banded minhash
            overhead = node_overhead(minhash)
            containment = node_containment(minhash)
            if containment and overhead <= max_overhead:
                add_node(node_id, minhash)
                return
        else:                             # empty minhash: do we keep?
            is_full, containment, overhead = var_in_bf_decide(node_id)

            # do we guesstimate this is good? ok => keep.
            if containment and overhead <= max_overhead:
                add_node(node_id, None)
                return

            # is there any chance a good hash is below this? if not, exit.
            if not is_full and containment == 0:
                return

        children_ids = dag[node_id]

        # leaf node. good varhash? keep.
        if not children_ids:
            if var_in_bf(node_id):
                add_node(node_id, minhash)
            return

        # recurse into children
        for child_id in children_ids:
            add_to_frontier2(child_id)

    add_to_frontier2(top_node_id)

    print('visited:', len(seen_nodes))

    return frontier, num_leaves, num_empty, frontier_minhash


def main(args=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('query_sig', help='query minhash')
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('overhead', help='\% of overhead', type=float)
    p.add_argument('--no-empty', action='store_true')
    p.add_argument('--purgatory', action='store_true')
    p.add_argument('-o', '--output', default=None)
    p.add_argument('--fullstats', action='store_true')
    p.add_argument('--checkfrontier', action='store_true')
    p.add_argument('-k', '--ksize', default=31, type=int,
                        help='k-mer size (default: 31)')

    args = p.parse_args(args)

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')

    # load catlas DAG
    top_node_id, dag, dag_up, dag_levels, cdbg_to_catlas = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    # load minhash DB
    db_path = get_minhashdb_name(args.catlas_prefix, args.ksize, 0, 0, 42)
    if not db_path:
        print('** ERROR, minhash DB does not exist for {}'.format(args.ksize),
              file=sys.stderr)
        sys.exit(-1)
    minhash_db = search_utils.MinhashSqlDB(db_path)

    # load query MinHash
    query_sig = load_query_signature(args.query_sig, ksize=args.ksize,
                                     select_moltype='DNA')
    print('loaded query sig {}'.format(query_sig.name()))

    frontier, num_leaves, num_empty, frontier_mh = frontier_search(query_sig, top_node_id, dag, minhash_db, args.overhead, not args.no_empty, args.purgatory)

    top_mh = load_minhash(top_node_id, minhash_db)

    scaled = max(top_mh.scaled, query_sig.minhash.scaled)
    query_mh = query_sig.minhash.downsample_scaled(scaled)
    top_mh = top_mh.downsample_scaled(scaled)

    print("Root containment: {}".format(query_mh.contained_by(top_mh)))
    print("Root similarity: {}".format(query_mh.similarity(top_mh)))

    print("Containment of frontier: {}".format(query_mh.contained_by(frontier_mh)))
    print("Similarity of frontier: {}".format(query_mh.similarity(frontier_mh)))
    print("Size of frontier: {} of {} ({:.3}%)".format(len(frontier), len(dag), 100 * len(frontier) / len(dag)))
    print("Overhead of frontier: {}".format(compute_overhead(frontier_mh, query_mh)))
    print("Number of leaves in the frontier: {}".format(num_leaves))
    print("Number of empty catlas nodes in the frontier: {}".format(num_empty))
    print("")

    if args.fullstats:
        shadow = find_shadow(frontier, dag)

        print("Size of the cDBG shadow: {}".format(len(shadow)))

    if args.checkfrontier:
        test_mh = top_mh.copy_and_clear()
        for node_id in frontier:
            mh = load_minhash(node_id, minhash_db)
            if mh:
                mh = mh.downsample_scaled(test_mh.scaled)
                test_mh.merge(mh)

        assert test_mh.similarity(frontier_mh) == 1.0

    query_size = len(query_sig.minhash.get_mins())
    query_bp = query_size * query_sig.minhash.scaled
    print("Size of query minhash: {} (est {:2.1e} bp)".\
              format(query_size, query_bp))
    minhash_size = len(frontier_mh.get_mins())
    minhash_bp = minhash_size * frontier_mh.scaled
    print("Size of frontier minhash: {} (est {:2.1e} bp); ratio {:.2f}".\
              format(minhash_size, minhash_bp, minhash_bp / query_bp))

    if args.output:
        print('saving frontier minhash as sourmash signature, into {}'.format(args.output))
        with open(args.output, 'w') as fp:
            sig = signature.SourmashSignature('', frontier_mh,
                                              name='frontier o={:1.2f} {}'.format(args.overhead, str(args.output)))
            sourmash_lib.signature.save_signatures([sig], fp)

    log(args.catlas_prefix, sys.argv)


if __name__ == '__main__':
    # import cProfile
    # cProfile.run('main()', 'search_stats')

    main()
