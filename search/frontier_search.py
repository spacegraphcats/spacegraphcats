#! /usr/bin/env python
import argparse
import os
import sys
from copy import copy
import heapq

import sourmash_lib
from sourmash_lib import MinHash, signature
from sourmash_lib.sourmash_args import load_query_signature
from typing import Dict, List, Set, Union, Tuple

from .memoize import memoize
from . import search_utils
from .search_utils import load_dag
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


def frontier_search_exact(top_node_id, dag, node_kmer_sizes, node_query_kmers, max_overhead):

    # nodes and minhases that are in the frontier
    frontier = []  # type: List[int]
    frontier_minhash = None

    # nodes and minhashes for leaves with large overhead
    seen_nodes = set()  # type: Set[int]

    num_leaves = 0
    num_empty = 0

    def node_containment(node_id):
        query_kmers = node_query_kmers.get(node_id, 0)
        total_kmers = node_kmer_sizes[node_id]
        containment = query_kmers / total_kmers
        assert query_kmers <= total_kmers

        return containment

    def node_overhead(node_id):
        return 1 - node_containment(node_id)

    def add_node(node_id: int):
        frontier.append(node_id)

    n_truncated = 0

    # This visits every node for which containment != 1.0, and adds
    # all children with any kind of containment.
    def add_to_frontier4(node_id):
        nonlocal n_truncated

        if node_id in seen_nodes:
            return

        seen_nodes.add(node_id)

        containment = node_containment(node_id)
        if containment == 0:
            n_truncated += 1
            return

        if containment >= 1.0 - max_overhead:
            add_node(node_id)
            n_truncated += 1
            return

        children_ids = dag[node_id]

        # leaf node: some containment? keep.
        if not children_ids and containment > 0.0:
            add_node(node_id)
            return

        # recurse into children to get more resolution
        for child_id in children_ids:
            add_to_frontier4(child_id)

    add_to_frontier4(top_node_id)

    frontier = set(frontier)

    print('search visited {} catlas nodes.'.format(len(seen_nodes)))
    print('search truncated at {}'.format(n_truncated))

    return frontier, num_leaves, num_empty, frontier_minhash


def frontier_search(top_node_id,
                    dag,
                    node_kmer_sizes,
                    node_query_kmers,
                    max_overhead,
                    min_containment):
    # the initial frontier is just the root
    root_match = node_query_kmers.get(top_node_id, 0)
    root_total = node_kmer_sizes[top_node_id]
    root_precision = root_match/root_total

    # since python implements heaps as min-heaps but we want to prioritize the
    # largest overheads, we store the precision (matches/total), since that is
    # 1 - overhead
    frontier = [(root_precision, root_match, root_total, top_node_id)]
    f_match = root_match
    f_total = root_total

    # keep refining until we either get the frontier overhead low enough, the
    # containment gets too low, or we run out of vertices
    while frontier and\
            1.0 - f_match/f_total > max_overhead and\
            f_match/root_match > min_containment:
        # remove worst overhead from the frontier
        _, curr_match, curr_total, curr = heapq.heappop(frontier)
        # update the quality metrics to reflect removal
        f_total = f_total - curr_total
        f_match = f_match - curr_match
        # descend into children
        for child in dag[curr]:
            child_match = node_query_kmers.get(child, 0)
            # if the child has no match, we don't want to include it
            if child_match == 0:
                continue

            # compute the other metrics
            child_total = node_kmer_sizes[child]
            child_precision = child_match/child_total
            f_total += child_total
            f_match += child_match

            # add to frontier
            heapq.heappush(frontier,
                           (child_precision, child_match, child_total, child))

    if len(frontier) > 0:
        frontier = list(zip(*frontier))[3]
    return frontier, 0, 0, None


def main(args):
    pass


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('--overhead', help='\% of overhead', type=float,
                   default=0.0)
    p.add_argument('output')
    p.add_argument('--query', help='query sequences', nargs='+')
    p.add_argument('--no-empty', action='store_true')
    p.add_argument('-k', '--ksize', default=31, type=int,
                   help='k-mer size (default: 31)')
    p.add_argument('--scaled', default=1000, type=float)
    p.add_argument('-v', '--verbose', action='store_true')
