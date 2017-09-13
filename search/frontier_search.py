#! /usr/bin/env python
import argparse
import os
import sys
import leveldb
from copy import copy

import sourmash_lib
from sourmash_lib import MinHash, signature
from sourmash_lib.sourmash_args import load_query_signature
from typing import Dict, List, Set, Union, Tuple

from .memoize import memoize
from .search_catlas_with_minhash import load_dag, load_minhash
from spacegraphcats.logging import log


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


def frontier_search(query_sig, top_node_id: int, dag, minhash_db: Union[str, leveldb.LevelDB], max_overhead: float, include_empty = True, use_purgatory = False):
    """
        include_empty: Whether to include nodes with no minhases in the frontier
        use_purgatory: put leaf nodes with large overhead into a purgatory and consider them after we have processes the whole frontier
    """
    # expand the frontier where the child nodes together have more than x% overhead

    # load the leveldb unless we get a path
    if not isinstance(minhash_db, leveldb.LevelDB):
        minhash_db = leveldb.LevelDB(minhash_db)

    # nodes and minhases that are in the frontier
    frontier = []  # type: List[int]
    frontier_minhash = None

    # nodes and minhashes for leaves with large overhead
    purgatory = []  # type: List[Tuple[float, int, MinHash]]

    seen_nodes = set()  # type: Set[int]

    num_leaves = 0
    num_empty = 0

    @memoize
    def load_and_downsample_minhash(node_id: int) -> MinHash:
        minhash = load_minhash(node_id, minhash_db)
        if minhash is None:
            return None
        return minhash.downsample_max_hash(query_sig.minhash)

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
        frontier.append(node_id)
        if minhash:
            if frontier_minhash:
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
            raise Exception("No containment")

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

    add_to_frontier(top_node_id)

    # now check whether the nodes in the purgatory are still necessary
    if len(purgatory):
        purgatory.sort()
        required_query_minhashes = query_sig.minhash.subtract_mins(frontier_minhash)
        for _, node_id, node_mh in purgatory:
            before = len(required_query_minhashes)
            required_query_minhashes -= set(node_mh.get_mins())
            after = len(required_query_minhashes)
            if before > after:
                frontier_minhash.merge(node_mh)
                frontier.append(node_id)
                num_leaves += 1

            if len(required_query_minhashes) == 0:
                # early termination, the full query is covered
                break

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
    p.add_argument('-k', '--ksize', default=None, type=int,
                        help='k-mer size (default: 31)')

    args = p.parse_args(args)

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')

    # load catlas DAG
    top_node_id, dag, dag_levels = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    db_path = os.path.join(args.catlas_prefix, 'minhashes.db')
    minhash_db = leveldb.LevelDB(db_path)

    # load query MinHash
    query_sig = load_query_signature(args.query_sig, select_ksize=args.ksize,
                                     select_moltype='DNA')
    print('loaded query sig {}'.format(query_sig.name()))

    frontier, num_leaves, num_empty, frontier_mh = frontier_search(query_sig, top_node_id, dag, minhash_db, args.overhead, not args.no_empty, args.purgatory)

    top_mh = load_minhash(top_node_id, minhash_db)
    query_mh = query_sig.minhash.downsample_max_hash(top_mh)
    top_mh = top_mh.downsample_max_hash(query_sig.minhash)
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
