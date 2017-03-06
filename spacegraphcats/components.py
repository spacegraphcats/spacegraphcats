from spacegraphcats.Eppstein import UnionFind
from collections import defaultdict

def component_index(graph):
    """
    Returns a union-find data structure associating each vertex 
    with an id (another vertex) the connected component it is contained in.
    """
    comps = UnionFind()
    for v in graph:
        N = graph.in_neighbors(v, 1).union([v])
        comps.union(*N)
    return comps

def components(graph):
    comps = component_index(graph)
    res = defaultdict(set)
    for v in graph:
        res[comps[v]].add(v)
    return res.values()

def num_components(graph) -> int:
    """ Compute the number of components """
    comps = component_index(graph)
    ids = set([comps[v] for v in comps])
    return len(ids)
