import sys, itertools, operator
from collections import defaultdict as defaultdict
from Eppstein import priorityDictionary, UnionFind
import hashlib

class VertexDict(dict):
    @classmethod
    def from_vxt(cls, file, param_parser=None):
        from parser import _parse_line
        if param_parser == None:
            param_parser = lambda p: p

        res = cls()
        for line in file:
            line = _parse_line(line)
            u, params = line[0], line[1:]
            res[int(u)] = param_parser(params)
        return res

    @classmethod
    def from_mxt(cls, file):
        from minhash import MinHash
        from parser import _parse_line
        res = cls()
        for line in file:
            u, hashes = _parse_line(line)
            hashes = list(map(int,hashes.split()))
            res[int(u)] = MinHash.from_list(hashes)
        return res


    def write_vxt(self, file, param_writer=None):
        if param_writer == None:
            param_writer = lambda p: p

        for u, p in self.items():
            file.write('{},{}\n'.format(u, param_writer(p)))


class EdgeStream:
    @staticmethod 
    def from_ext(file):
        from parser import _parse_line
        for line in file:
            u, v = list(map(int, _parse_line(line)))
            yield (u,v)

class EdgeSet(set):
    @classmethod 
    def from_ext(cls, file):
        from parser import _parse_line
        res = cls()
        for u, v in EdgeStream.from_ext(file):
            res.add((u, v))

        return res

    def write_ext(self, file):
        for u, v in self:
            file.write('{},{}\n'.format(u,v))


class Graph:
    def __init__(self):
        self.adj = defaultdict(set)
        self.nodes = set()

    @staticmethod
    def from_gxt(file):
        from parser import parse
        res = Graph()
        node_attr = defaultdict(dict)
        edge_attr = defaultdict(dict)
        def add_node(v, size, prop_names, props):
            res.add_node(v)
            node_attr[v]['size'] = size
            for key,value in zip(prop_names,props):
                node_attr[v][key] = value

        def add_edge(u, v, prop_names, props):
            res.add_edge(u,v)
            for key,value in zip(prop_names,props):
                edge_attr[(u,v)][key] = value

        parse(file, add_node, add_edge)
        return res, node_attr, edge_attr

    @staticmethod
    def on(nodes):
        res = Graph()
        res.nodes = set(nodes)
        return res

    def __contains__(self,u):
        return u in self.nodes

    def __iter__(self):
        return iter(self.nodes)

    def __len__(self):
        return len(self.nodes)

    def get_max_id(self):
        return max(self.nodes)

    def edges(self):
        for u in self:
            for v in self.adj[u]:
                if u <= v:
                    yield (u,v)

    def num_edges(self):
        return sum(1 for _ in self.edges())

    def add_node(self,u):
        self.nodes.add(u)
        return self

    def add_edge(self,u,v):
        self.nodes.add(u)
        self.nodes.add(v)
        self.adj[u].add(v)
        self.adj[v].add(u)
        return self

    def remove_edge(self,u,v):
        self.adj[u].discard(v)
        self.adj[v].discard(u)
        return self

    def remove_node(self, u):
        self.nodes.remove(u)
        for v in self.adj[u]:
            self.adj[v].remove(u)
        del self.adj[u]
        return self        

    def remove_loops(self):
        for v in self:
            self.remove_edge(v,v)
        return self            

    def has_loops(self):
        for v in self:
            if self.adjacent(v,v):
                return True
        return False

    def adjacent(self,u,v):
        return v in self.adj[u]

    def neighbours(self,u):
        return frozenset(self.adj[u])

    def neighbours_set(self, centers):
        return self.neighbours_set_closed(centers) - centers

    def neighbours_set_closed(self, centers):
        res = set()
        for v in centers:
            res = res | self.neighbours(v)
        return res

    def rneighbours(self,u,r):
        res = set([u])
        for _ in range(r):
            res |= self.neighbours_set_closed(res)
        return res

    def degree(self,u):
        return len(self.adj[u])

    def degree_sequence(self):
        return [ self.degree(u) for u in self.nodes]

    def degree_dist(self):
        res = defaultdict(int)
        for u in self.nodes:
            res[self.degree(u)] += 1
        return res

    def calc_average_degree(self):
        num_edges = len( [e for e in self.edges()] )
        return 2.0*num_edges / len(self.nodes)

    def calc_degeneracy(self):
        degbuckets = defaultdict(set)
        degrees = {}
        degen = 0
        mindeg = n = len(self)
        removed = set()
        for v in self:
            deg = self.degree(v)
            mindeg = min(mindeg, deg)
            degbuckets[deg].add(v)
            degrees[v] = deg

        while len(removed) < n:
            while len(degbuckets[mindeg]) == 0:
                mindeg += 1
            v = next(iter(degbuckets[mindeg]))
            degbuckets[mindeg].remove(v)
            removed.add(v)
            degen = max(degen, mindeg)

            for w in self.neighbours(v):
                if w in removed:
                    continue
                deg = degrees[w]
                degrees[w] -= 1
                degbuckets[deg].remove(w)
                degbuckets[deg-1].add(w)
                mindeg = min(mindeg, deg-1)

        return degen

    def copy(self):
        res = Graph()
        for v in self:
            res.add_node(v)
        for u,v in self.edges():
            res.add_edge(u,v)
        return res

    def subgraph(self, vertices):
        vertices = frozenset(vertices)
        res = Graph.on(vertices)
        for u in vertices:
            N = self.neighbours(u) & vertices
            for v in N:
                res.add_edge(u,v)
        return res

    # Returns graph with node indices from 0 to n-1 and
    # a mapping dictionary to recover the original ids
    def normalize(self):
        res = Graph()
        backmapping = dict(enumerate(self))
        mapping = dict( (k,v) for (v,k) in backmapping.items() )

        for u in self:
            res.nodes.add( mapping[u] )
        for u,v in self.edges():
            res.add_edge( mapping[u], mapping[v] )

        return res, backmapping

    # For memoization
    def __str__(self):
        return graph_hash(self)

    def is_connected(self):
        if len(self) <= 1:
            return True
        comp = set()
        comp.add(next(iter(self)))
        exp = self.neighbours_set(comp)
        while len(exp) > 0:
            comp.update( exp )
            exp = self.neighbours_set(comp)

        return len(self) == len(comp)

    def component_index(self):
        """ Returns a union-find data structure associating each vertex 
            with an id (another vertex) the connected component it is contained in """
        comps = UnionFind()
        for v in self:
            N = self.neighbours(v) | set([v])
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

class TFGraph:
    def __init__(self, nodes):
        self.nodes = nodes
        self.inarcs = defaultdict(dict)
        self.inarcs_weight = defaultdict(lambda : defaultdict(set))

    def __contains__(self,u):
        return u in self.nodes

    def __len__(self):
        return len(self.nodes)

    def __iter__(self):
        return iter(self.nodes)

    def add_node(self,u):
        self.nodes.add(u)
        return self

    def add_arc(self, u, v, weight):
        assert u != v
        self.inarcs[v][u] = weight
        self.inarcs_weight[v][weight].add(u)
        return self

    def add_arcs(self, arcs, weight):
        for u, v in arcs:
            self.add_arc(u, v, weight)
        return self

    def remove_arc(self,u,v):
        if u not in self.inarcs[v]:
            return
        weight = self.inarcs[v][u]
        del self.inarcs[v][u]
        self.inarcs_weight[v][weight].remove(u)
        assert not self.adjacent(u,v)
        return self

    def arcs(self, weight=None):
        if weight == None:
            for u in self:
                for v, w in self.in_neighbours(u):
                    yield (v, u, w)
        else:
            for u in self:
                for v, w in self.in_neighbours(u):
                    if w == weight:
                        yield (v, u)

    def num_arcs(self):
        return sum( 1 for _ in self.arcs() )

    # For compatibility reasons (coloring etc.)
    def edges(self):
        for u in self:
            for v,_ in self.in_neighbours(u):
                yield (v,u)

    # Return the arc weight of uv or None if it does not exit
    def weight(self,u,v):
        return self.inarcs[v][u] if u in self.inarcs[v] else None

    # Returns whether the arc uv exists
    def adjacent(self,u,v):
        return u in self.inarcs[v]

    def in_neighbours(self,u):
        inbs = self.inarcs[u]
        for v in inbs:
            yield v, inbs[v]

    def in_neighbours_weight(self,u,weight):
        inbs = self.inarcs_weight[u][weight]
        for v in inbs:
            yield v

    def in_degree(self,u):
        return len(self.inarcs[u])

    def in_degree_weight(self,u):
        res = {}
        for weight in self.inarcs_weight[u]:
            res[weight] = len(self.inarcs_weight[u][weight])
        return res

    def in_degree_sequence(self):
        return [ self.in_degree(u) for u in self.nodes]

    def in_degree_dist(self):
        res = defaultdict(int)
        for u in self.nodes:
            res[self.in_degree(u)] += 1
        return res

    def degree_dist(self):
        degrees = defaultdict(int)
        for u in self.nodes:
            for w in self.inarcs[u]:
                degrees[u] += 1
                degrees[w] += 1
        res = defaultdict(int)
        for u in self.nodes:
            res[degrees[u]] += 1
        return res

    def distance(self,u,v):
        distance = float('inf')
        Nu = {}
        for x,d in self.in_neighbours(u):
            Nu[x] = d
            if x == v:
                distance = d
        for x,d in self.in_neighbours(v):
            if x == u:
                distance = min(distance,d)
            elif x in Nu:
                distance = min(distance,d+Nu[x])
        return distance

    def copy(self):
        res = TFGraph(self.nodes)
        for u,v,d in self.arcs():
            res.add_arc(u,v,d)
        return res

    # Returns an undirected copy
    def undirected(self):
        res = Graph()
        for v in self.nodes:
            res.nodes.add(v) # For degree-0 nodes!
        for u,v,_ in self.arcs():
            res.add_edge(u,v)
        return res

    # Returns transitive triples (x,u,weightsum)
    def trans_trips(self, u):
        inbs = frozenset(self.inarcs[u].items())
        for (y, wy) in inbs:
            assert y != u
            for x, wx in self.in_neighbours(y):
                if not self.adjacent(x, u):
                    yield (x, u, wx+wy)

    # Returns transitive triples (x,y,weightsum) whose weightsum
    # is precisely 'weight'
    def trans_trips_weight(self,u,weight):
        for wy in range(1, weight):
            wx = weight-wy
            inbs = frozenset(self.inarcs_weight[u][wy])
            for y in inbs:
                assert y != u
                for x in self.in_neighbours_weight(y, wx):
                    if not self.adjacent(x, u) and x != u:
                        yield (x, u, weight)


    # Returns fraternal triples (x,y,weightsum)
    def frat_trips(self,u):
        inbs = frozenset(self.inarcs[u].items())
        for (x, wx), (y, wy) in itertools.combinations(inbs, 2):
            assert x != u and y != u
            if not (self.adjacent(x, y) or self.adjacent(y, x)):
                yield (x, y, wx+wy)

    # Returns fraternal triples (x,y,weightsum) whose weightsum
    # is precisely 'weight'
    def frat_trips_weight(self, u, weight):
        ''' Since we draw all vertices from the same in-neighbourhood,
        the considered sets have usually different weights and
        are thus disjoint. Only for even weights do we have to consider
        the slightly annoying case of (weight/2, weight/2).
        '''
        wh = (weight+1) // 2
        for wx in range(1,wh):
            wy = weight-wx
            for x in self.in_neighbours_weight(u,wx):
                assert x != u
                for y in self.in_neighbours_weight(u,wy):
                    assert y != u
                    if not (self.adjacent(x,y) or self.adjacent(y,x)):
                        yield (x, y, weight)

        if weight % 2 == 0:
            wy = wx = weight // 2
            inbs = frozenset(self.inarcs_weight[u][wh])
            for x, y in itertools.combinations(inbs, 2):
                assert x != u and y != u
                if not (self.adjacent(x,y) or self.adjacent(y,x)):
                    yield (x, y, weight)

    # For memoization
    def __str__(self):
        return graph_hash(self)


def write_gxt(file, g, node_attrs=None, edge_attrs=None):
    from parser import Writer

    node_keys = set()
    if node_attrs:
        for u in g:
            node_keys |= set(node_attrs[u].keys())

    size_key = None
    if 'size' in node_keys:
        size_key = 'size'
    elif 'segmentsize' in node_keys:
        size_key = 'segmentsize'
    node_keys.discard(size_key)
    node_keys = sorted(list(node_keys))

    edge_keys = set()
    if edge_attrs:
        for u,v in g.edges():
            edge_keys |= set(edge_attrs[(u,v)])

    edge_keys = sorted(list(edge_keys))

    writer = Writer(file, node_keys, edge_keys)

    if size_key:
        for u in g:
            attributes = [node_attrs[u][key] if key in node_attrs[u] else None for key in node_keys]
            writer.add_vertex(u,node_attrs[u][size_key],attributes)
    else:
        for u in g:
            attributes = [node_attrs[u][key] if key in node_attrs[u] else None for key in node_keys]
            writer.add_vertex(u,0,attributes)

    for u,v in g.edges():
            attributes = [edge_attrs[(u,v)][key] if key in edge_attrs[(u,v)] else None for key in edge_keys]
            writer.add_edge(u,v,attributes)
    writer.done()

def graph_hash(g):
    if isinstance(g, TFGraph):
        deg = defaultdict(int)
        for v in g:
            indegs = g.in_degree_weight(v)
            for w in indegs:
                deg[ (indegs[w],w) ] += 1
        m = hashlib.md5()
        for d in sorted(deg.keys()):
            m.update(str(deg[d]**(d[0]*d[1])).encode())
        return m.hexdigest()
    elif isinstance(g, Graph):
        deg = defaultdict(int)
        for v in g:
            deg[g.degree(v)] += 1
        m = hashlib.md5()
        for d in sorted(deg.keys()):
            m.update(str(deg[d]**d).encode())
        return m.hexdigest()
