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
