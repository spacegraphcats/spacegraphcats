import itertools
from spacegraphcats.Eppstein import UnionFind
from collections import defaultdict

class Graph(object):
    """
    Graph data structure.  Undirected edges are represented by bidirectional
     directed arcs
    """
    def __init__(self, num_nodes, radius=1):
        self.num_nodes = num_nodes
        self.radius = radius
        # lookup the weight of the arc by its two endpoints
        self.inarcs = [dict() for _ in range(self.num_nodes)]
        # lookup the neighbors of a vertex with a given weight
        self.inarcs_by_weight = [[set() for i in range(self.num_nodes)] 
                                  for j in range(self.radius)]

    def __iter__(self):
        return iter(range(self.num_nodes))

    def __len__(self):
        return self.num_nodes

    def add_arc(self, u, v, weight=1):
        self.inarcs[v][u] = weight
        # use weight-1 for 0 indexing
        self.inarcs_by_weight[weight-1][v].add(u)
        return self

    def remove_arc(self, u, v):
        """
        remove arc uv from the graph
        Precondition/assumption:  arc uv has a unique weight
        """
        if u not in self.inarcs[v]:
            return
        # find the weight of the arc
        weight = self.inarcs[v][u]
        del self.inarcs[v][u]
        self.inarcs_by_weight[weight-1][v].remove(u)
        return self

    def adjacent(self, u, v):
        """
        check if uv or vu is an arc
        """
        return u in self.inarcs[v] or v in self.inarcs[u]

    def arcs(self,weight=None):
        """
        return all the arcs in the graph.  restrict to a given weight when provided
        """
        if weight:
            return self.inarcs_by_weight[weight]
        else:
            return self.inarcs

    def in_neighbors(self, v, weight=None):
        if weight is None:
            return self.inarcs[v].items()
        else:
            return self.inarcs_by_weight[weight-1][v]

    def in_degree(self, v, weight=None):
        return len(self.in_neighbors(v, weight))

    def num_arcs(self):
        return sum(1 for _ in self.arcs())

    def transitive_pairs(self, u, w):
        """
        Returns transitive pairs (x,u) whose weight sum is exactly w
        """
        # loop over all weights that sum to w
        for wy in range(1, w):
            wx = w-wy
            for y in self.in_neighbors(u,wy):
                #assert y != u
                for x in self.in_neighbors(y, wx):
                    # need to double check that the x,u edges isn't there already
                    if x != u and not self.adjacent(x, u):
                        yield x, u


    def fraternal_pairs(self,u,w):
        """
        Returns fraternal pairs of in_neighbors of u (x,y) whose weight sum is
        exactly w.  (x,y) are fraternal if they are not adjacent.
        """
        # Since we draw all vertices from the same in-neighborhood,
        # the considered sets have usually different weights and
        # are thus disjoint. Only for even weights do we have to consider
        # the slightly annoying case of (w/2, w/2).
        wh = (w+1) // 2
        for wx in range(1,wh):
            wy = w-wx
            for x in self.in_neighbors(u,wx):
                #assert x != u
                for y in self.in_neighbors(u,wy):
                    #assert y != u
                    if not (self.adjacent(x,y) or self.adjacent(y,x)):
                        yield x, y
        # slightly annoying case
        if w % 2 == 0:
            inbs = self.in_neighbors(u,wh)
            for x, y in itertools.combinations(inbs, 2):
                assert x != u and y != u
                if not (self.adjacent(x,y) or self.adjacent(y,x)):
                    yield x, y 

    def component_index(self):
        """
        Returns a union-find data structure associating each vertex 
        with an id (another vertex) the connected component it is contained in.
        """
        comps = UnionFind()
        for v in self:
            N = self.in_neighbors(v, 1) | set([v])
            comps.union(*N)
        return comps

    def components(self):
        comps = self.component_index()
        res = defaultdict(set)
        for v in self:
            res[comps[v]].add(v)
        return res.values()

    def num_components(self):
        comps = self.component_index()
        ids = set([comps[v] for v in comps])
        return len(ids)

class DictGraph(Graph):
    def __init__(self, nodes=None, r=1):
        if nodes is None:
            self.nodes = set()
        else:
            self.nodes = nodes
        self.radius = r
        self.inarcs = defaultdict(dict)
        self.inarcs_by_weight = [defaultdict(set) for _ in range(self.radius)]

    def __len__(self):
        return len(self.nodes)

    def __iter__(self):
        return iter(self.nodes)

    def add_node(self,u):
        self.nodes.add(u)
        return self
