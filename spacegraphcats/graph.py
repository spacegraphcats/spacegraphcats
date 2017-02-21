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
        # lookup the neighbors of a vertex with a given weight
        self.inarcs_by_weight = [[set() for i in range(self.num_nodes)] 
                                  for j in range(self.radius)]

    def __iter__(self):
        return iter(range(self.num_nodes))

    def __len__(self):
        return self.num_nodes

    def add_arc(self, u, v, weight=1):
        # use weight-1 for 0 indexing
        self.inarcs_by_weight[weight-1][v].add(u)
        return self

    def remove_arc(self, u, v):
        """
        remove arc uv from the graph
        Precondition/assumption:  arc uv has a unique weight
        """
        for arc_list in self.inarcs_by_weight:
            if u in arc_list[v]:
                arc_list[v].remove(u)
                return self
        return self

    def adjacent(self, u, v):
        """
        check if uv or vu is an arc
        """
        for arc_list in self.inarcs_by_weight:
            if u in arc_list[v] or v in arc_list[u]:
                return True
        return False

    def arcs(self,weight=None):
        """
        return all the arcs in the graph.  restrict to a given weight when provided
        """
        if weight is None:
            return [(y,x,w+1) for w,arc_list in enumerate(self.inarcs_by_weight) 
                for x,N in enumerate(arc_list)
                for y in N]
        else:
            return [(y,x) for x,N in enumerate(self.inarcs_by_weight[weight-1]) 
                for y in N]

    def in_neighbors(self, v, weight=None):
        if weight is None:
            #return self.inarcs[v].items()
            return [(u,w+1) for w,arc_list in enumerate(self.inarcs_by_weight) 
                for u in arc_list[v]]
        else:
            return self.inarcs_by_weight[weight-1][v]

    def in_degree(self, v, weight=None):
        if weight is None:
            return sum(len(x[v]) for x in self.inarcs_by_weight)
        else:
            return len(self.inarcs_by_weight[weight-1][v])

    def num_arcs(self, weight=None):
        return sum(self.in_degree(v, weight) for v in self)

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
        self.inarcs_by_weight = [defaultdict(set) for _ in range(self.radius)]

    def __len__(self):
        return len(self.nodes)

    def __iter__(self):
        return iter(self.nodes)

    def add_node(self,u):
        self.nodes.add(u)

    def arcs(self,weight=None):
        """
        return all the arcs in the graph.  restrict to a given weight when provided
        """
        if weight:
            return [(x,y) for x,N in self.inarcs_by_weight[weight-1].items() 
                for y in N]
        else:
            return [(x,y,w+1) for w,arc_list in enumerate(self.inarcs_by_weight) 
                for x,N in arc_list.items()
                for y in N]