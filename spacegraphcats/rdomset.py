#!/usr/bin/env python
from collections import defaultdict, Counter
from operator import itemgetter
import itertools
import sys, random, os
import argparse
from Eppstein import priorityDictionary
from graph import Graph, TFGraph, write_gxt

""" Dtf functions """

def dtf(g, r):
    auggraph = ldo(g)
    num_arcs = auggraph.num_arcs()

    changed = True
    d = 1
    print "Augmenting",
    while changed and d <= r:
        print d,
        dtf_step(auggraph, d+1)

        curr_arcs = auggraph.num_arcs() # This costs a bit so we store it
        changed = num_arcs < curr_arcs
        num_arcs = curr_arcs
        d += 1
    print ""
    return auggraph

def dtf_step(g, dist):
    fratGraph = Graph()
    newTrans = {}

    for v in g:
        for x,y,_ in g.trans_trips_weight(v, dist):
            newTrans[(x,y)] = dist
        for x,y,_ in g.frat_trips_weight(v, dist):
            fratGraph.add_edge(x,y)

    for (s, t) in newTrans:
        g.add_arc(s, t, dist)
        fratGraph.remove_edge(s,t)

    fratDigraph = ldo(fratGraph)

    for s,t,_ in fratDigraph.arcs():
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
    for i in xrange(0, len(g)):
        d = 0
        while len(buckets[d]) == 0:
            d += 1
        v = iter(buckets[d]).next()
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
    for _ in xrange(2):
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
    for _ in xrange(2):
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
        for r in xrange(d+1):
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


def calc_domset_graph(domset, dominators, d):
    h = Graph.on(domset)

    for v in dominators:
        if v in domset:
            continue
        # Collect closest pairs
        sorted_doms = filter(lambda s: len(s) > 0, [dominators[v][r] for r in xrange(d+1)])
        if len(sorted_doms[0]) >= 2:
            for x,y in itertools.combinations(sorted_doms[0],2):
                h.add_edge(x,y)
        else:
            assert len(sorted_doms[0]) == 1
            x = iter(sorted_doms[0]).next()
            if len(sorted_doms) >= 2:
                for y in sorted_doms[1]:
                    h.add_edge(x,y)

    return h

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

def vc_compression(g,dominators,d):
    confgraph = Graph()
    forced = set()
    for v in g:
        candidates = [x for r in xrange(d+1) for x in dominators[v][r] ]
        if len(candidates) == 1:
            confgraph.add_node(candidates[0])
            forced.add(candidates[0])
            continue
        elif len(candidates) == 2:
            confgraph.add_edge(*candidates)
            continue

        fixed = [] # Candidate already in the graph
        for c in reversed(candidates):
            if c in confgraph:
                fixed.append(c)
                if len(fixed) == 2:
                    break
        u,w = None,None
        if len(fixed) == 0:
            u,w = random.sample(candidates,2)
        elif len(fixed) == 1:
            u = fixed[0]
            w = random.sample(candidates,1)[0]
        else:
            u,w = fixed
        confgraph.add_edge(u,w)

    # Calculate vertex cover of confgraph
    cover = set(forced)
    edges = [(u,v,g.degree(u),g.degree(v)) for u,v in confgraph.edges()]
    edges.sort(key=lambda (u,v,du,dv): du+dv, reverse=False)
    for u,v,_,_ in edges:
        if u not in cover and v not in cover:
            cover.add(u)
            cover.add(v)
    domset = cover

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
        for r in xrange(0,d+1):
            candidates -= dominators[v][r]
        if len(candidates) > 0:
            print v, candidates, dominators[v]
            return False
    return True

def test_dominators_strict(g, domset, dominators, d):
    for v in g:
        dN = set([v])
        if dominators[v][0] != (dN & domset):
            return False

        for r in xrange(1,d+1):
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
    print "Computed {}-domset of size {} in graph with {} vertices".format(d,len(domset),len(g))

    dominators = calc_dominators(augg, domset, d)

    # Compute domset graph
    if domgraph:
        print "Computing domgraph"
        h = calc_domset_graph(domset, dominators, d)
        write_gxt(domgraph, h)

    # Test domset
    if test:
        print "Testing domset:", test_domset(g,domset,d)

    # Annotate domset
    labelds = 'ds{}'.format(d)
    for v in domset:
        node_attrs[v][labelds] = 1

    write_gxt(outputgxt, g, node_attrs, edge_attrs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='gml input file', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('output', help='gml output file', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdin)
    parser.add_argument('--domgraph', help='gml output file', nargs='?', type=argparse.FileType('w'),
                        default=None)
    parser.add_argument("d", help="domset distance", type=int )
    parser.add_argument('--test', action='store_true')
    args = parser.parse_args()

    main(args.input, args.output, args.d, args.test, args.domgraph)

