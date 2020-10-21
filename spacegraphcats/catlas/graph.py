"""Graph data structure."""
import itertools
from collections import defaultdict

from typing import List, Set, Iterable, Any, Sized

from sortedcontainers import SortedSet


class Graph(Iterable, Sized):
    """
    Graph data structure.

    Undirected edges are represented by bidirectional directed arcs
    """

    def __init__(self, num_nodes: int, radius: int = 1) -> None:
        """Create a new graph."""
        self.num_nodes = num_nodes
        self.radius = radius
        # lookup the neighbors of a vertex with a given weight
        self.inarcs_by_weight = [
            [SortedSet() for i in range(self.num_nodes)] for j in range(self.radius)
        ]  # type: Any

    def __iter__(self):
        """Iteration support."""
        return iter(range(self.num_nodes))

    def __len__(self):
        """len() support."""
        return self.num_nodes

    def __contains__(self, v):
        """Support for `if v in graph`."""
        return v < self.num_nodes

    def add_arc(self, u: int, v: int, weight: int = 1):
        """Add the arc u,v at a given weight."""
        # use weight-1 for 0 indexing
        self.inarcs_by_weight[weight - 1][v].add(u)
        return self

    def remove_arc(self, u: int, v: int):
        """
        Remove arc uv from the graph.

        Precondition/assumption:  arc uv has a unique weight
        """
        for arc_list in self.inarcs_by_weight:
            if u in arc_list[v]:
                arc_list[v].remove(u)
                return self
        return self

    def adjacent(self, u: int, v: int):
        """Check if uv or vu is an arc."""
        for arc_list in self.inarcs_by_weight:
            if u in arc_list[v] or v in arc_list[u]:
                return True
        return False

    def arcs(self, weight: int = None):
        """
        Return an iterator over all the arcs in the graph.

        Restrict to a given weight when provided
        """
        if weight is None:
            for w, arc_list in enumerate(self.inarcs_by_weight):
                for x, N in enumerate(arc_list):
                    for y in N:
                        yield y, x, w + 1
        else:
            for x, N in enumerate(self.inarcs_by_weight[weight - 1]):
                for y in N:
                    yield y, x

    def in_neighbors(self, v: int, weight: int = None):
        """Return the inneighbors at a vertex."""
        if weight is None:
            # return self.inarcs[v].items()
            return [
                (u, w + 1)
                for w, arc_list in enumerate(self.inarcs_by_weight)
                for u in arc_list[v]
            ]
        else:
            return self.inarcs_by_weight[weight - 1][v]

    def in_degree(self, v: int, weight: int = None):
        """Return the number of in neighbors."""
        if weight is None:
            return sum(len(x[v]) for x in self.inarcs_by_weight)
        else:
            return len(self.inarcs_by_weight[weight - 1][v])

    def num_arcs(self, weight: int = None):
        """Return the total number of arcs in this graph."""
        return sum(self.in_degree(v, weight) for v in self)

    def transitive_pairs(self, u: int, w: int):
        """Return transitive pairs (x,u) whose weight sum is exactly w."""
        # loop over all weights that sum to w
        for wy in range(1, w):
            wx = w - wy
            for y in self.in_neighbors(u, wy):
                # assert y != u
                for x in self.in_neighbors(y, wx):
                    # need to double check that the x,u edges isn't there
                    # already
                    if x != u and not self.adjacent(x, u):
                        yield x, u

    def fraternal_pairs(self, u: int, w: int):
        """
        Return fraternal pairs of in_neighbors of u (x,y) whose weight sum is
        exactly w.  (x,y) are fraternal if they are not adjacent.
        """
        # It is sufficient to consider pairs whose weights are 1,w-1.  When w
        # > 2, 1 and w-1 are distinct and their neighborhoods are thus
        # disjoint.  When w == 2, we need to ignore pairs of the form (x,x)
        if w == 2:
            inbs = self.in_neighbors(u, 1)
            for x, y in itertools.combinations(inbs, 2):
                assert x != u and y != u
                if not self.adjacent(y, x):
                    yield x, y
        else:
            for x in self.in_neighbors(u, 1):
                # assert x != u
                for y in self.in_neighbors(u, w - 1):
                    # assert y != u
                    if not self.adjacent(y, x):
                        yield x, y


class DictGraph(Graph):
    """Graph that supports nonconsecutive vertex ids."""

    def __init__(self, nodes: Set[int] = None, r: int = 1) -> None:
        """Make a new graph."""
        if nodes is None:
            self.nodes = SortedSet()  # type: Set[int]
        else:
            self.nodes = nodes
        self.radius = r
        self.inarcs_by_weight = [defaultdict(SortedSet) for _ in range(self.radius)]

    def __len__(self):
        """len() support."""
        return len(self.nodes)

    def __iter__(self):
        """Iteration support."""
        return iter(self.nodes)

    def __contains__(self, v):
        """Support for `if v in graph`."""
        return v in self.nodes

    def add_node(self, u: int):
        """Add a new node."""
        self.nodes.add(u)

    def arcs(self, weight: int = None):
        """
        Return all the arcs in the graph.

        restrict to a given weight when provided
        """
        if weight:
            return [
                (x, y) for x in self.nodes for y in self.inarcs_by_weight[weight - 1][x]
            ]
        else:
            return [
                (x, y, w + 1)
                for w, arc_list in enumerate(self.inarcs_by_weight)
                for x in self.nodes
                for y in arc_list[x]
            ]

    def remove_isolates(self) -> List[int]:
        """
        Remove all isolated vertices and return a list of vertices removed.

        Precondition:  the graph is bidirectional
        """
        isolates = [v for v in self if self.in_degree(v, 1) == 0]
        for v in isolates:
            self.nodes.remove(v)
            for N in self.inarcs_by_weight:
                if v in N:
                    del N[v]
        return isolates
