#!/usr/bin/env python3
from __future__ import print_function
import itertools, sys, random, os, argparse, glob
from collections import defaultdict, Counter
from operator import itemgetter
from os import path
from Eppstein import priorityDictionary
from graph import Graph, TFGraph, EdgeSet, EdgeStream, VertexDict, write_gxt


class Domination:
    def __init__(self, domset, domgraph, assignment, radius):
        self.domset = domset
        self.domgraph = domgraph
        self.assignment = VertexDict(assignment)
        self.radius = radius

    def write(self, projectpath, projectname):
        fname = path.join(projectpath, projectname+".{name}.{radius}.{extension}")

        with open(fname.format(name="domgraph",extension="gxt",radius=self.radius), 'w') as f:
            write_gxt(f, self.domgraph)

        with open(fname.format(name="assignment",extension="vxt",radius=self.radius), 'w') as f:
            self.assignment.write_vxt(f, param_writer=lambda s: ' '.join(map(str,s)))

    @staticmethod
    def read(projectpath, projectname, radius):
        fname = path.join(projectpath, projectname+".{name}.{radius}.{extension}")

        with open(fname.format(name="domgraph",extension="gxt",radius=radius), 'r') as f:
            domgraph, _, _ = Graph.from_gxt(f)

        with open(fname.format(name="assignment",extension="vxt",radius=radius), 'r') as f:
            assignment = VertexDict.from_vxt(f, lambda s: list(map(int,s[0].split())))

        domset = set([v for v in domgraph])

        return Domination(domset, domgraph, assignment, radius)


class LazyDomination:
    def __init__(self, project):
        self.radius = project.radius
        self.projectpath = project.path
        self.projectname = project.name
        self.graph = project.graph

    def compute(self):
        try:
            res = Domination.read(self.projectpath, self.projectname, self.radius)
            print("Loaded {}-domination from project folder".format(self.radius))
        except IOError:
            augg = self._compute_augg()
            domset = better_dvorak_reidl(augg, self.radius)
            dominators = calc_dominators(augg, domset, self.radius)
            domgraph, domset, dominators, assignment = calc_domination_graph(self.graph, augg, domset, dominators, self.radius)

            res = Domination(domset, domgraph, assignment, self.radius)
            res.write(self.projectpath, self.projectname)

        return res

    def _compute_augg(self):
        """
            Computes dtf augmentations or loads them from the project folder.
        """
        augname = path.join(self.projectpath,self.projectname+".aug.{}.ext")

        augs = {}
        for f in glob.glob(augname.format("[0-9]*")):
            d = int(f.split(".")[-2])
            augs[d] = f

        if 0 in augs:
            auggraph = TFGraph(self.graph)
            with open(augname.format("0"), 'r') as f:
                auggraph.add_arcs(EdgeStream.from_ext(f), 1)
        else:
            auggraph = ldo(self.graph)
            with open(augname.format("0"), 'w') as f:
                EdgeSet(auggraph.arcs(weight=1)).write_ext(f)

        num_arcs = auggraph.num_arcs()
        changed = True
        d = 1
        print("Augmenting", end=" ")
        sys.stdout.flush()
        while changed and d <= self.radius:
            if d in augs:
                print("({})".format(d), end=" ")
                sys.stdout.flush()
                with open(augname.format(d), 'r') as f:
                    auggraph.add_arcs(EdgeStream.from_ext(f), d+1)
            else:
                print(d, end=" ")
                sys.stdout.flush()
                dtf_step(auggraph, d+1)
                with open(augname.format(d), 'w') as f:
                    EdgeSet(auggraph.arcs(weight=d+1)).write_ext(f)            

            curr_arcs = auggraph.num_arcs() # This costs a bit so we store it
            changed = num_arcs < curr_arcs
            num_arcs = curr_arcs
            d += 1
        print("")
        return auggraph



""" Dtf functions """

def dtf(g, r):
    auggraph = ldo(g)
    num_arcs = auggraph.num_arcs()

    changed = True
    d = 1
    print("Augmenting", end=" ", flush=True)
    while changed and d <= r:
        print(d, end=" ", flush=True)
        dtf_step(auggraph, d+1)

        curr_arcs = auggraph.num_arcs() # This costs a bit so we store it
        changed = num_arcs < curr_arcs
        num_arcs = curr_arcs
        d += 1
    print("", flush=True)
    return auggraph

def dtf_step(g, dist):
    fratGraph = Graph()
    newTrans = {}

    for v in g:
        for x, y, _ in g.trans_trips_weight(v, dist):
            assert x != y
            newTrans[(x, y)] = dist
        for x, y, _ in g.frat_trips_weight(v, dist):
            assert x != y
            fratGraph.add_edge(x, y)

    for (s, t) in newTrans:
        assert s != t
        g.add_arc(s, t, dist)
        fratGraph.remove_edge(s,t)

    fratDigraph = ldo(fratGraph)

    for s, t, _ in fratDigraph.arcs():
        assert s != t
        g.add_arc(s,t,dist)

def ldo(g, weight=None):
    res = TFGraph(g.nodes)

    if weight == None:
        weight = defaultdict(int)

    degdict = {}
    buckets = defaultdict(set)
    for v in g:
        d = g.degree(v) + weight[v]
        degdict[v] = d
        buckets[d].add(v)

    seen = set()
    for i in range(0, len(g)):
        d = 0
        while len(buckets[d]) == 0:
            d += 1
        v = next(iter(buckets[d]))
        buckets[d].remove(v)

        for u in g.neighbours(v):
            if u in seen:
                continue
            d = degdict[u]
            buckets[d].remove(u)
            buckets[d-1].add(u)
            degdict[u] -= 1
            # Orient edges towards v
            res.add_arc(u,v,1)

        seen.add(v)
    return res

""" Dvorak's algorithm with dtf-augs instead of weakly linked cols """

def calc_distance(augg, X):
    dist = defaultdict(lambda: float('inf'))

    for v in X:
        dist[v] = 0

    # We need two passes
    for _ in range(2):
        for v in augg:
            for u, r in augg.in_neighbours(v):
                dist[v] = min(dist[v],dist[u]+r)
            for u, r in augg.in_neighbours(v):
                dist[u] = min(dist[u],dist[v]+r)

    return dist

def calc_dominators(augg, domset, d):
    dominators = defaultdict(lambda: defaultdict(set))

    for v in domset:
        dominators[v][0].add(v)

    # We need two passes
    for _ in range(2):
        for v in augg:
            # Pull dominators
            for u, r in augg.in_neighbours(v):
                for r2 in dominators[u].keys():
                    domdist = r+r2
                    if domdist <= d:
                        dominators[v][domdist] |= dominators[u][r2]
            # Push dominators
            for u, r in augg.in_neighbours(v):
                for r2 in dominators[v].keys():
                    domdist = r+r2
                    if domdist <= d:
                        dominators[u][domdist] |= dominators[v][r2]

    for v in dominators:
        cumulated = set()
        for r in range(d+1):
            dominators[v][r] -= cumulated
            cumulated |= dominators[v][r]

    return dominators

def calc_dominator_sizes(g, domset, d):
    res = {}
    for v in domset:
        res[v] = len(g.rneighbours(v,d))
    return res

def calc_dominated(augg, domset, d):
    dist = calc_distance(augg, domset)
    return filter(lambda v: dist[v] <= d, dist)


def _calc_domset_graph(domset, dominators, d):
    """ 
        Builds up a 'domination graph' by assigning each vertex to 
        its closest dominators. These dominators will be connected in the
        final graph.
    """
    h = Graph.on(domset)

    assignment = defaultdict(set)
    for v in dominators:
        # Collect closest pairs
        sorted_doms = list(filter(lambda s: len(s) > 0, [dominators[v][r] for r in range(d+1)]))
        if len(sorted_doms[0]) >= 2:
            for x in sorted_doms[0]:
                assignment[v].add(x)

            for x,y in itertools.combinations(sorted_doms[0],2):
                h.add_edge(x,y)
        else:
            assert len(sorted_doms[0]) == 1
            x = next(iter(sorted_doms[0]))
            assignment[v].add(x)
            if len(sorted_doms) >= 2:
                for y in sorted_doms[1]:
                    assignment[v].add(y)
                    h.add_edge(x,y)

    return h, assignment

def calc_domination_graph(g, augg, domset, dominators, d):
    h, assignment = _calc_domset_graph(domset, dominators, d)
    hcomps = h.component_index()

    # print("  First domgraph has {} components".format(h.num_components()))

    compmap = {}
    for u in g:
        candidates = [hcomps[v] for v in assignment[u]]
        assert len(candidates) > 0
        assert candidates.count(candidates[0]) == len(candidates) # Make sure all the comps. agree
        compmap[u] = candidates[0]

    conn_graph = Graph()
    for u,v in g.edges():
        if u in domset or v in domset:
            continue
        # Check whether u,v 'connects' to components of h
        if compmap[u] != compmap[v]:
            conn_graph.add_edge(u,v)
 
    # Compute vertex cover for conn_graph
    vc = set()
    for u in conn_graph: # Add neighbours of deg-1 vertices
        if conn_graph.degree(u) == 1 and u not in vc:
            vc.add(next(iter(conn_graph.neighbours(u))))

    # We expect this to be mostly disjoint edges, so the following
    # heuristic works well.
    for u,v in conn_graph.edges():
        if u not in vc and v not in vc:
            if conn_graph.degree(u) > conn_graph.degree(v):
                vc.add(u)
            else:
                vc.add(v)

    newdomset = domset | vc
    if len(newdomset) != len(domset):
        # print("  Recomputing domgraph")
        dominators = calc_dominators(augg, newdomset, d)
        h, assignment = _calc_domset_graph(newdomset, dominators, d)

    return h, newdomset, dominators, assignment

def better_dvorak_reidl(augg,d):
    domset = set()
    infinity = float('inf')
    domdistance = defaultdict(lambda: infinity)

    #  low indegree first (reverse=False)
    vprops = [(v,augg.in_degree(v)) for v in augg]
    vprops.sort(key=itemgetter(1),reverse=False)
    order = map(itemgetter(0),vprops)

    for v in order:
        if domdistance[v] <= d:
            continue

        for u,r in augg.in_neighbours(v):
            domdistance[v] = min(domdistance[v],r+domdistance[u])

        if domdistance[v] <= d:
            continue

        domset.add(v)
        domdistance[v] = 0

        for u,r in augg.in_neighbours(v):
            domdistance[u] = min(domdistance[u],r+domdistance[v])

    return domset

def test_scattered(g, sc, d):
    sc = set(sc)
    for v in sc:
        dN = g.rneighbours(v,d)
        if len(dN & sc) > 1: # v is in dN but nothing else should
            return False
    return True

def test_domset(g, ds, d):
    dominated = set()
    for v in ds:
        dominated |= g.rneighbours(v,d)
    return len(dominated) == len(g)

def test_dominators(g, domset, dominators, d):
    for v in g:
        dN = g.rneighbours(v,d)
        candidates = dN & domset
        for r in range(0,d+1):
            candidates -= dominators[v][r]
        if len(candidates) > 0:
            print(v, candidates, dominators[v])
            return False
    return True

def test_dominators_strict(g, domset, dominators, d):
    for v in g:
        dN = set([v])
        if dominators[v][0] != (dN & domset):
            return False

        for r in range(1,d+1):
            shell = g.neighbours_set(dN) - dN
            if dominators[v][r] != (shell & domset):
                return False
            dN |= shell
    return True

def rdomset(g, d):
    augg = dtf(g, d)
    domset = better_dvorak_reidl(augg, d)

    return domset, augg

def main(inputgxt, outputgxt, d, test, domgraph):
    g, node_attrs, edge_attrs = Graph.from_gxt(inputgxt)
    g.remove_loops()

    domset, augg = rdomset(g,d)
    print("Computed {}-domset of size {} in graph with {} vertices".format(d,len(domset),len(g)))

    dominators = calc_dominators(augg, domset, d)

    # Compute domset graph
    if domgraph:
        print("Computing domgraph")
        h = _calc_domset_graph(domset, dominators, d)
        write_gxt(domgraph, h)

    # Test domset
    if test:
        print("Testing domset:", test_domset(g,domset,d))

    # Annotate domset
    labelds = 'ds{}'.format(d)
    for v in domset:
        node_attrs[v][labelds] = 1

    write_gxt(outputgxt, g, node_attrs, edge_attrs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='gxt input file', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('output', help='gxt output file', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdin)
    parser.add_argument('--domgraph', help='gxt output file', nargs='?', type=argparse.FileType('w'),
                        default=None)
    parser.add_argument("d", help="domset distance", type=int )
    parser.add_argument('--test', action='store_true')
    args = parser.parse_args()

    main(args.input, args.output, args.d, args.test, args.domgraph)

