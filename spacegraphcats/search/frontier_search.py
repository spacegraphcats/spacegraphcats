#! /usr/bin/env python
import heapq
import time
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Union


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


def frontier_search_exact(catlas,
                          node_query_kmers,
                          max_overhead):

    # nodes and minhases that are in the frontier
    frontier = []  # type: List[int]
    frontier_minhash = None

    # nodes and minhashes for leaves with large overhead
    seen_nodes = set()  # type: Set[int]

    num_leaves = 0
    num_empty = 0

    def node_containment(node_id):
        query_kmers = node_query_kmers.get(node_id, 0)
        total_kmers = catlas.index_sizes[node_id]
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

        children_ids = catlas.children[node_id]

        # leaf node: some containment? keep.
        if not children_ids and containment > 0.0:
            add_node(node_id)
            return

        # recurse into children to get more resolution
        for child_id in children_ids:
            add_to_frontier4(child_id)

    add_to_frontier4(catlas.root)

    frontier = set(frontier)

    print('search visited {} catlas nodes.'.format(len(seen_nodes)))
    print('search truncated at {}'.format(n_truncated))

    return frontier, num_leaves, num_empty, frontier_minhash


def frontier_search(catlas,
                    node_query_kmers,
                    max_overhead,
                    min_containment):
    # the initial frontier is just the root
    root_match = node_query_kmers.get(catlas.root, 0)
    root_total = catlas.index_sizes[catlas.root]
    root_precision = root_match/root_total

    # since python implements heaps as min-heaps but we want to prioritize the
    # largest overheads, we store the precision (matches/total), since that is
    # 1 - overhead
    frontier = [(root_precision, root_match, root_total, catlas.root)]
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
        for child in catlas.children[curr]:
            child_match = node_query_kmers.get(child, 0)
            # if the child has no match, we don't want to include it
            if child_match == 0:
                continue

            # compute the other metrics
            child_total = catlas.index_sizes[child]
            child_precision = child_match/child_total
            f_total += child_total
            f_match += child_match

            # add to frontier
            heapq.heappush(frontier,
                           (child_precision, child_match, child_total, child))

    if len(frontier) > 0:
        frontier = list(zip(*frontier))[3]
    return frontier, 0, 0, None


def collect_frontier(catlas,
                     catlas_match_counts,
                     max_overhead,
                     min_containment,
                     verbose=False):
    start = time.time()

    # gather results into total_frontier
    total_frontier = defaultdict(set)

    # do queries!
    try:
        frontier, num_leaves, num_empty, frontier_mh = \
          frontier_search(catlas,
                          catlas_match_counts,
                          max_overhead,
                          min_containment)
    except NoContainment:
        print('** WARNING: no containment!?')
        frontier = []

    # record which seed (always 0, here) contributed to which node
    for node in frontier:
        total_frontier[node].add(0)

    end = time.time()
    print('catlas query time: {:.1f}s'.format(end-start))

    return total_frontier


def collect_frontier_exact(catlas,
                           catlas_match_counts,
                           overhead=0.0,
                           verbose=False):
    start = time.time()

    # gather results into total_frontier
    total_frontier = defaultdict(set)

    # do queries!
    try:
        frontier, num_leaves, num_empty, frontier_mh = \
          frontier_search_exact(catlas,
                                catlas_match_counts,
                                overhead)

    except NoContainment:
        print('** WARNING: no containment!?')
        frontier = []

    # record which seed (always 0, here) contributed to which node
    for node in frontier:
        total_frontier[node].add(0)

    end = time.time()
    print('catlas query time: {:.1f}s'.format(end-start))

    return total_frontier
