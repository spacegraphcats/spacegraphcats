"""Compute the components of a graph."""
from collections import defaultdict
from typing import Set, ValuesView, Dict

from .Eppstein import UnionFind
from .graph import Graph


def component_index(graph: Graph) -> UnionFind:
    """
    Return a union-find data structure associating each vertex with an id
     (another vertex) the connected component it is contained in.
    """
    comps = UnionFind()
    for v in graph:
        N = graph.in_neighbors(v, 1).union([v])
        comps.union(*N)
    return comps


def components(graph: Graph) -> ValuesView[Set[int]]:
    """Return the components."""
    comps = component_index(graph)
    res = defaultdict(set)  # type: Dict[int, Set[int]]
    for v in graph:
        res[comps[v]].add(v)
    return res.values()


def num_components(graph: Graph) -> int:
    """Compute the number of components."""
    comps = component_index(graph)
    ids = set([comps[v] for v in comps])
    return len(ids)
