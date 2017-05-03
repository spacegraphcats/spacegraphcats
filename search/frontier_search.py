#! /usr/bin/env python
import os
import sys
import argparse
from sourmash_lib import MinHash
from collections import defaultdict
from spacegraphcats.catlas import CAtlas
import time
import sourmash_lib
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import search_minhashes, SigLeaf
from sourmash_lib import signature
import screed
import pickle
from typing import List

from .search_catlas_with_minhash import load_dag, load_minhash


def find_shadow(nodes: List[int], dag) -> List[int]:
    shadow = []

    def add_to_shadow(node_id):
        children_ids = dag[node_id]

        if len(children_ids) == 0:
            shadow.append(node_id)
        else:
            for child in children_ids:
                add_to_shadow(child)
    
    for node in nodes:
        add_to_shadow(node)

    return shadow


def compute_overhead(node_minhash, query_minhash) -> float:
    """ Compute the relative overhead of minhashes. """
    node_length = len(node_minhash.get_mins())
    return (node_length - node_minhash.count_common(query_minhash)) / node_length


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query_sig', help='query minhash')
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('overhead', help='\% of overhead', type=float)

    args = p.parse_args()

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, basename + '.catlas')

    # load catlas DAG
    top_node_id, dag, dag_levels = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    # load query MinHash
    query_sig = sourmash_lib.signature.load_signatures(args.query_sig)
    query_sig = list(query_sig)[0]
    print('loaded query sig {}'.format(query_sig.name()))

    minhash_dir = os.path.join(args.catlas_prefix, basename + '.minhashes')

    # expand the frontier where the child nodes together have more than x% overhead

    frontier = []
    seen_nodes = set()

    def add_to_frontier(node_id: int):
        """
        Add a node or its children to the frontier.
        """
        if (node_id in seen_nodes):
            # we have added this node or its children to the frontier
            return
        
        seen_nodes.add(node_id)

        children_ids = dag[node_id]

        if len(children_ids) == 0:
            # leaf
            frontier.append(node_id)
            return

        # check whether the node has more than x% overhead
        minhash = load_minhash(node_id, minhash_dir)

        query_mh = query_sig.minhash.downsample_max_hash(minhash)
        against_mh = minhash.downsample_max_hash(query_sig.minhash)

        containment = query_mh.containment(against_mh)

        if containment == 0:
            # ignore nodes that have no containment
            # sanity check, should not get here
            raise Exception("No containment")

        overhead = compute_overhead(against_mh, query_mh)

        print("{} Overhead {}".format(node_id, overhead))
 
        if overhead > args.overhead:
            # greedily find the minimum number of children such that they together contain everything the node contains

            overheads = []

            for child in children_ids:
                child_mh = load_minhash(child, minhash_dir)

                # ignore children without any containment
                if query_mh.containment(child_mh) == 0:
                    continue

                child_oh = compute_overhead(child_mh, query_mh)
                overheads.append((child_oh, child, child_mh))

            overheads.sort()

            print(overheads)

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
            print("low overhead")
            frontier.append(node_id)

    add_to_frontier(top_node_id)

    
    union = load_minhash(frontier.pop(), minhash_dir)
    for node in frontier:
        mh = load_minhash(node, minhash_dir)
        union.merge(mh)

    print(union.similarity(query_sig.minhash))

    print("Size of frontier: {} of {} ({:.2}%)".format(len(frontier), len(dag), 100.0 * len(frontier) / len(dag)))

    shadow = find_shadow(frontier, dag)

    print("Size shadow: {}".format(len(shadow)))

    sys.exit(0)


if __name__ == '__main__':
    main()
