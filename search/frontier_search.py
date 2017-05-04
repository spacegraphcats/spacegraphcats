#! /usr/bin/env python
import argparse
import os
import sys

import sourmash_lib
from sourmash_lib import MinHash, signature
from typing import Dict, List, Set

from .memoize import memoize
from .search_catlas_with_minhash import load_dag, load_minhash


def find_shadow(nodes: List[int], dag: Dict[int, List[int]]) -> Set[int]:
    shadow = set()  # type: Set[int]

    seen_nodes = set()  # type: Set[int]

    def add_to_shadow(node_id: int):
        if node_id in seen_nodes:
            return

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
    """ Compute the relative overhead of minhashes. """
    node_length = len(node_minhash.get_mins())
    return (node_length - node_minhash.count_common(query_minhash)) / node_length


def frontier_search(query_sig, top_node_id: int, dag, minhash_dir: str, max_overhead: float):
    # expand the frontier where the child nodes together have more than x% overhead

    query_mh = query_sig.minhash
    frontier = []
    frontier_minhash = None
    seen_nodes = set()  # type: Set[int]

    num_leaves = 0

    @memoize
    def load_and_downsample_minhash(node_id: int):
        nonlocal query_mh
        minhash = load_minhash(node_id, minhash_dir)
        if minhash:
            minhash = minhash.downsample_max_hash(query_mh)
            query_mh = query_mh.downsample_max_hash(minhash)
        return minhash

    @memoize
    def node_overhead(minhash):
        return compute_overhead(minhash, query_mh)

    @memoize
    def node_containment(minhash):
        return query_mh.containment(minhash)

    def add_node(node_id: int, minhash):
        nonlocal frontier_minhash
        frontier.append(node_id)
        if frontier_minhash:
            frontier_minhash.merge(minhash)
        else:
            frontier_minhash = minhash


    def add_to_frontier(node_id: int):
        """
        Add a node or its children to the frontier.
        """
        if node_id in seen_nodes:
            # we have added this node or its children to the frontier
            return
        else:
            seen_nodes.add(node_id)

        children_ids = dag[node_id]

        minhash = load_and_downsample_minhash(node_id)
            
        if len(children_ids) == 0:
            # leaf
            nonlocal num_leaves
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

            for child in children_ids:
                child_mh = load_and_downsample_minhash(child)
                if not child_mh:
                    continue

                # ignore children without any containment
                if node_containment(child_mh) == 0:
                    continue

                child_oh = node_overhead(child_mh)
                overheads.append((child_oh, child, child_mh))

            overheads.sort()

            first = overheads.pop()
            add_to_frontier(first[1])
            union = first[2]
            for child_tuple in overheads:
                # add children to frontier
                add_to_frontier(child_tuple[1])

                union.merge(child_tuple[2])
                if union.containment(query_mh) == containment:
                    # early termination, all children already cover the node so we can stop
                    return

        else:
            # low overhead node
            add_node(node_id, minhash)

    add_to_frontier(top_node_id)

    return frontier, num_leaves, frontier_minhash


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query_sig', help='query minhash')
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('overhead', help='\% of overhead', type=float)
    p.add_argument('-o', '--output', default=None)

    args = p.parse_args()

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, basename + '.catlas')

    # load catlas DAG
    top_node_id, dag, dag_levels = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    minhash_dir = os.path.join(args.catlas_prefix, basename + '.minhashes')

    # load query MinHash
    query_sig = sourmash_lib.signature.load_signatures(args.query_sig)
    query_sig = list(query_sig)[0]
    print('loaded query sig {}'.format(query_sig.name()))

    frontier, num_leaves, frontier_mh = frontier_search(query_sig, top_node_id, dag, minhash_dir, args.overhead)

    top_mh = load_minhash(top_node_id, minhash_dir)
    query_mh = query_sig.minhash.downsample_max_hash(top_mh)
    top_mh = top_mh.downsample_max_hash(query_sig.minhash)
    print("Root containment: {}".format(query_mh.containment(top_mh)))
    print("Root similarity: {}".format(query_mh.similarity(top_mh)))

    print("Containment of frontier: {}".format(query_mh.containment(frontier_mh)))
    print("Similarity of frontier: {}".format(query_mh.similarity(frontier_mh)))
    print("Size of frontier: {} of {} ({:.3}%)".format(len(frontier), len(dag), 100 * len(frontier) / len(dag)))
    print("Number of leaves in the frontier: {}".format(num_leaves))

    shadow = find_shadow(frontier, dag)

    print("Size of the shadow: {}".format(len(shadow)))

    if args.output:
        print('saving frontier minhash as sourmash signature, into {}'.format(args.output))
        with open(args.output, 'w') as fp:
            sig = signature.SourmashSignature('', frontier_mh,
                                              name='frontier')
            sourmash_lib.signature.save_signatures([sig], fp)

    sys.exit(0)


if __name__ == '__main__':
    # import cProfile
    # cProfile.run('main()', 'search_stats')

    main()
