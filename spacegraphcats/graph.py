import sys, itertools, operator
from collections import defaultdict as defaultdict
from Eppstein import priorityDictionary, UnionFind
import hashlib

class Graph(object):
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

    def add_edge(self,u,v):
        self.nodes.add(u)
        self.nodes.add(v)
        self.adj[u].add(v)
        self.adj[v].add(u)

    def remove_edge(self,u,v):
        self.adj[u].discard(v)
        self.adj[v].discard(u)

    def remove_node(self, u):
        self.nodes.remove(u)
        for v in self.adj[u]:
            self.adj[v].remove(u)
        del self.adj[u]

    def remove_loops(self):
        for v in self:
            self.remove_edge(v,v)

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
        for _ in xrange(r):
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
        return float(2*num_edges) / len(self.nodes)

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
            v = iter(degbuckets[mindeg]).next()
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
        res = Graph()
        selected = set(vertices)
        for v in self:
            if v in selected:
                res.add_node(v)

        for u,v in self.edges():
            if u in selected and v in selected:
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

    @staticmethod
    def generate(name):
        res = Graph()

        # Ki, Pi, Ci
        t = name[0]
        if t in ['K','P','C']:
            if t == 'K':
                parts = map( int, name[1:].split(',') )
                if len(parts) == 1:
                    size = parts[0]
                    for x,y in itertools.combinations(xrange(0,size),2):
                        res.add_edge(x,y)
                else:
                    print("Generating", name)
                    indices = []
                    size = 0
                    for s in parts:
                        indices.append( (size, size+s) )
                        size += s
                    for a,b in itertools.combinations(xrange(0,len(parts)),2):
                        for x in xrange( indices[a][0], indices[a][1] ):
                            for y in xrange( indices[b][0], indices[b][1] ):
                                res.add_edge(x,y)

            elif t == 'P' or t == 'C':
                size = int(name[1:])
                for x in xrange(0,size-1):
                    res.add_edge(x,x+1)
                if t == 'C':
                    res.add_edge(size-1,0)

            return res

        raise Exception("Unknown name "+name)

    # For memoization
    def __str__(self):
        return graph_hash(self)

    def is_connected(self):
        if len(self) <= 1:
            return True
        comp = set( [iter(self).next()] )
        exp = self.neighbours_set(comp)
        while len(exp) > 0:
            comp.update( exp )
            exp = self.neighbours_set(comp)

        return len(self) == len(comp)

    def get_components(self, vertices=None):
        if vertices == None:
            vertices = set(self.nodes)
        else:
            vertices = set(vertices)

        while len(vertices) > 0:
            comp = set([ vertices.pop()])
            exp = self.neighbours_set(comp) & vertices
            while len(exp) > 0:
                comp.update( exp )
                exp = self.neighbours_set(comp) & vertices

            vertices = vertices - comp
            yield comp

def graph_from_pickle(lst):
    cls = lst.pop(0) # load class object
    g = Graph()

    while len(lst) > 0:
        u = lst.pop(0)
        g.add_edge(u,lst.pop(0))

    return g

class TFGraph(object):
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

    def add_arc(self,u,v,weight):
        self.inarcs[v][u] = weight
        self.inarcs_weight[v][weight].add(u)

    def remove_arc(self,u,v):
        if u not in self.inarcs[v]:
            return
        weight = self.inarcs[v][u]
        del self.inarcs[v][u]
        self.inarcs_weight[v][weight].remove(u)
        assert not self.adjacent(u,v)

    def arcs(self):
        for u in self:
            for v,weight in self.in_neighbours(u):
                yield (v,u,weight)

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
    def trans_trips(self,u):
        inbs = frozenset(self.inarcs[u].items())
        for (y,wy) in inbs:
            for x,wx in self.in_neighbours(y):
                if not self.adjacent(x,u):
                    yield (x,u,wx+wy)

    # Returns transitive triples (x,y,weightsum) whose weightsum
    # is precisely 'weight'
    def trans_trips_weight(self,u,weight):
        for wy in xrange(1,weight):
            wx = weight-wy
            inbs = frozenset(self.inarcs_weight[u][wy])
            for y in inbs:
                for x in self.in_neighbours_weight(y,wx):
                    if not self.adjacent(x,u):
                        yield (x,u,weight)


    # Returns fraternal triples (x,y,weightsum)
    def frat_trips(self,u):
        inbs = frozenset(self.inarcs[u].items())
        for (x,wx),(y,wy) in itertools.combinations(inbs, 2):
            if not (self.adjacent(x,y) or self.adjacent(y,x)):
                yield (x,y,wx+wy)

    # Returns fraternal triples (x,y,weightsum) whose weightsum
    # is precisely 'weight'
    def frat_trips_weight(self,u,weight):
        ''' Since we draw all vertices from the same in-neighbourhood,
        the considered sets have usually different weights and
        are thus disjoint. Only for even weights do we have to consider
        the slightly annoying case of (weight/2, weight/2).
        '''
        wh = (weight+1)/2
        for wx in xrange(1,wh):
            wy = weight-wx
            for x in self.in_neighbours_weight(u,wx):
                for y in self.in_neighbours_weight(u,wy):
                    if not (self.adjacent(x,y) or self.adjacent(y,x)):
                        yield (x,y,weight)

        if weight % 2 == 0:
            wy = wx = weight/2
            inbs = frozenset(self.inarcs_weight[u][wh])
            for x,y in itertools.combinations(inbs, 2):
                if not (self.adjacent(x,y) or self.adjacent(y,x)):
                    yield (x,y,weight)

    # For memoization
    def __str__(self):
        return graph_hash(self)

class AGraph(object):
    def __init__(self):
        self.adj = defaultdict(set)
        self.nodes = set()
        self.degs = defaultdict(int)
        self.deglists = defaultdict(set)
        self.hidden = set()
        self.m = 0

    @staticmethod
    def from_networkx(nxg):
        res = AGraph()
        for u in nxg:
            res.add_node(u)
        for u,v in nxg.edges():
            res.add_edge(u,v)
        return res

    def add_node(self,u):
        if u in self.nodes:
            return
        self.nodes.add(u)
        self.degs[u] = 0
        self.deglists[0].add(u)

    def hide_node(self,u):
        assert u in self.nodes
        assert u not in self.hidden
        affected = set()
        for v in self.adj[u]: # Also decrease degree of hidden nodes!
            olddeg = self.degs[v]
            self.degs[v] -= 1
            self.m -= 0 if v in self.hidden else 1
            self.deglists[olddeg].remove(v)
            self.deglists[olddeg-1].add(v)
            affected.add(v)
        self.hidden.add(u)
        return affected - self.hidden

    def hide_nodes(self,nodes):
        for v in nodes:
            if v not in self.hidden:
                self.hide_node(v)

    def show_node(self,u):
        for v in self.adj[u]: # Also increase degree of hidden nodes!
            olddeg = self.degs[v]
            self.degs[v] += 1
            self.m += 0 if v in self.hidden else 1
            self.deglists[olddeg].remove(v)
            self.deglists[olddeg+1].add(v)
        self.hidden.remove(u)

    def show_nodes(self,nodes):
        for v in nodes:
            if v in self.hidden:
                self.show_node(v)

    def degree(self, u):
        assert u not in self.hidden
        return self.degs[u]

    def neighbours(self, u):
        return self.adj[u] - self.hidden

    def _add_adj(self,u,v):
        if v not in self.adj[u]:
            self.adj[u].add(v)
            # Update u's degree
            olddeg = self.degs[u]
            self.degs[u] = olddeg+1
            self.deglists[olddeg].remove(u)
            self.deglists[olddeg+1].add(u)
            if u < v:
                self.m += 1

    def add_edge(self,u,v):
        assert u not in self.hidden, v not in self.hidden
        self.add_node(u)
        self.add_node(v)
        self._add_adj(u,v)
        self._add_adj(v,u)

    def num_edges(self):
        return self.m

    def get_by_degree(self,d):
        return self.deglists[d] - self.hidden

    def num_components(self):
        comps = UnionFind()
        for v in self:
            N = self.neighbours(v) | set([v])
            comps.union(*N)
        ids = set([comps[v] for v in comps])
        return len(ids)

    def __len__(self):
        return len(self.nodes) - len(self.hidden)

    def __iter__(self):
        return iter(self.nodes - self.hidden)


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
            m.update(str(deg[d]**(d[0]*d[1])))
        return m.hexdigest()
    elif isinstance(g, Graph):
        deg = defaultdict(int)
        for v in g:
            deg[g.degree(v)] += 1
        m = hashlib.md5()
        for d in sorted(deg.keys()):
            m.update(str(deg[d]**d))
        return m.hexdigest()