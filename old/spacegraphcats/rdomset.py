#!/usr/bin/env python3
from __future__ import print_function

import itertools, sys, random, os, argparse, glob
from collections import defaultdict, Counter
from operator import itemgetter
from os import path

from spacegraphcats.Eppstein import priorityDictionary
from spacegraphcats.graph import Graph, DictGraph, TFGraph, EdgeSet, EdgeStream, VertexDict, write_gxt
from spacegraphcats.graph_parser import IdentityHash


class Domination:
    def __init__(self, domset, domgraph, assignment, radius):
        self.domset = domset
        self.domgraph = domgraph
        self.assignment = VertexDict(assignment)
        self.radius = radius

    def write(self, projectpath, projectname, id_map=None):
        fname = path.join(projectpath, projectname+".{name}.{radius}.{extension}")

        if id_map is None:
            id_map = IdentityHash()

        with open(fname.format(name="domgraph",extension="gxt",radius=self.radius), 'w') as f:
            write_gxt(f, self.domgraph, id_map=id_map)

        with open(fname.format(name="assignment",extension="vxt",radius=self.radius), 'w') as f:
            self.assignment.write_vxt(f, param_writer=lambda s: ' '.join(map(lambda x:str(id_map[x]),s)), id_map=id_map)

    @staticmethod
    def read(projectpath, projectname, radius, id_map = None):
        fname = path.join(projectpath, projectname+".{name}.{radius}.{extension}")

        if id_map is None:
            id_map = IdentityHash()

        with open(fname.format(name="domgraph",extension="gxt",radius=radius), 'r') as f:
            domgraph, _, _, _ = DictGraph.from_gxt(f)

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
        self.id_map = project.id_map

    def compute(self):
        try:
            res = Domination.read(self.projectpath, self.projectname, self.radius, self.id_map)
            print("Loaded {}-domination from project folder".format(self.radius))
        except IOError:
            augg = self._compute_augg()
            domset = better_dvorak_reidl(augg, self.radius)
            dominators = calc_dominators(augg, domset, self.radius)
            domgraph, assignment = calc_domination_graph(self.graph, augg, domset, dominators, self.radius)

            res = Domination(domset, domgraph, assignment, self.radius)
            res.write(self.projectpath, self.projectname, self.id_map)

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
            with open(augname.format("0"), 'r') as f:
                self.graph.add_arcs(EdgeStream.from_ext(f, self.id_map), 1)
        else:
            ldo(self.graph)
            with open(augname.format("0"), 'w') as f:
                EdgeSet(self.graph.arcs(weight=1)).write_ext(f,self.id_map)

        num_arcs = auggraph.num_arcs()
        changed = True
        d = 1
        print("Augmenting", end=" ")
        sys.stdout.flush()
        while changed and d < self.radius:
            if d in augs:
                print("({})".format(d), end=" ")
                sys.stdout.flush()
                with open(augname.format(d), 'r') as f:
                    auggraph.add_arcs(EdgeStream.from_ext(f,self.id_map), d+1)
            else:
                print(radius, end=" ")
                sys.stdout.flush()
                dtf_step(auggraph, d+1)
                with open(augname.format(d), 'w') as f:
                    EdgeSet(auggraph.arcs(weight=d+1)).write_ext(f,self.id_map)            

            curr_arcs = auggraph.num_arcs() # This costs a bit so we store it
            changed = num_arcs < curr_arcs
            num_arcs = curr_arcs
            d += 1
        print("")
        return auggraph


def dtf(g, radius, comp=None):
    """ Computes the r-th dft-augmentation of g. """
    auggraph = ldo(g,comp=comp)
    num_arcs = auggraph.num_arcs()

    changed = True
    d = 1
    while changed and d < radius:
        dtf_step(auggraph, d+1, comp)

        # Small optimization: if no new arcs have been added we can stop.
        curr_arcs = auggraph.num_arcs() # This costs a bit so we store it
        changed = num_arcs < curr_arcs
        num_arcs = curr_arcs
        d += 1
    return auggraph

def dtf_step(augg, dist,comp=None):
    """ Computes the d-th dtf-augmentation from a dtf-graph augg
        (where d is provided by the argument dist). The input augg
        must be a (d-1)-th dtf-augmentation. See dtf() for usage. 
        This function adds arcs to augg. """
    fratGraph = DictTFGraph() # Records fraternal edge, must be oriented at the end
    newTrans = {} # Records transitive arcs

    # if a list of nodes in a component is supplied, we loop over that,
    # otherwise we loop over all nodes in the graph.
    if comp is None:
        nodes = augg
    else:
        nodes = comp

    for v in nodes:
        for x, y, _ in augg.trans_trips_weight(v, dist):
            #assert x != y
            newTrans[(x, y)] = dist
        for x, y, _ in augg.frat_trips_weight(v, dist):
            #assert x != y
            fratGraph.add_arc(x, y)
            fratGraph.add_arc(y, x)
    # Add transitive arcs to graph
    for (s, t) in newTrans:
        #assert s != t
        augg.add_arc(s, t, dist)
        fratGraph.remove_arc(s,t)
        fratGraph.remove_arc(t,s)

    # Orient fraternal edges and add them to the grah
    ldo(fratGraph)

    for s, t, _ in fratDigraph.arcs():
        #assert s != t
        augg.add_arc(s,t,dist)

def ldo(g, weight=None, comp=None):
    """ Computes a low-in-degree orientation of a tf graph g in place
        by iteratively removing a vertex of mimimum degree and orienting
        the edges towards it. """

    # if a list of nodes in a component is supplied, we loop over that,
    # otherwise we loop over all nodes in the graph.
    if comp is None:
        nodes = g
    else:
        nodes = comp

    if weight == None:
        weight = defaultdict(int)

    degdict = {}
    buckets = defaultdict(set)
    for v in nodes:
        d = g.in_degree(v) + weight[v]
        degdict[v] = d
        buckets[d].add(v)

    seen = set()
    for i in range(0, len(nodes)):
        d = 0
        while len(buckets[d]) == 0:
            d += 1
        v = buckets[d].pop()

        for u, _ in g.in_neighbours(v):
            if u in seen:
                continue
            d = degdict[u]
            buckets[d].remove(u)
            buckets[d-1].add(u)
            degdict[u] -= 1
            # Orient edges towards v by removing the edge from v to u
            g.remove_arc(v,u)

        seen.add(v)

def calc_distance(augg, X, comp=None):
    """ Computes the for every vertex in the graph underlying the
        augmentation augg its distance to the set X """
    dist = defaultdict(lambda: float('inf'))

    if comp is None:
        nodes = augg
    else:
        nodes = comp

    # Vertices in X have distance zero to X
    for v in X:
        dist[v] = 0

    # We need two passes in order to only every access the
    # in-neighbourhoods (which are guaranteed to be small).
    for _ in range(2):
        for v in nodes:
            # Pull distances from in-neighbourhood
            for u, r in augg.in_neighbours(v):
                dist[v] = min(dist[v],dist[u]+r)
            # Push distances to in-neighbourhood                
            for u, r in augg.in_neighbours(v):
                dist[u] = min(dist[u],dist[v]+r)

    return dist

def calc_dominators(augg, domset, radius,comp=None):
    """
    Computes for each vertex the subset of domset that dominates it at distance d
    Returns a double level dictionary that maps vertices in the original graph to
    integers 0 to d to sets of vertices that dominate at that distance
    """
    dominators = defaultdict(lambda: defaultdict(set))

    if comp is None:
        nodes = augg
    else:
        nodes = comp

    # Every vertex in domset is a zero-dominator of itself
    for v in domset:
        dominators[v][0].add(v)

    # We need two passes in order to only every access the
    # in-neighbourhoods (which are guaranteed to be small).
    for _ in range(2):
        for v in nodes:
            # Pull dominators from in-neighbourhood
            for u, r in augg.in_neighbours(v):
                for r2 in dominators[u].keys():
                    domdist = r+r2
                    if domdist <= radius:
                        dominators[v][domdist] |= dominators[u][r2]
            # Push dominators to in-neighbourhood
            for u, r in augg.in_neighbours(v):
                for r2 in dominators[v].keys():
                    domdist = r+r2
                    if domdist <= radius:
                        dominators[u][domdist] |= dominators[v][r2]

    # Clean up: vertices might appear at multiple distances as dominators,
    # we only want to store them for the _minimum_ distance at which they
    # act as a dominator.
    for v in dominators:
        cumulated = set()
        for r in range(radius+1):
            dominators[v][r] -= cumulated
            cumulated |= dominators[v][r]

    return dominators

def calc_dominator_sizes(g, domset, radius):
    res = {}
    for v in domset:
        res[v] = len(g.rneighbours(v,radius))
    return res

def calc_dominated(augg, domset, radius):
    """ Compute the union of d-neighbourhoods around 
        the vertex set domset, that is, all vertices d-domainted
        by it. """
    dist = calc_distance(augg, domset)
    return filter(lambda v: dist[v] <= radius, dist)

def _calc_domset_graph(domset, dominators, radius):
    """ 
        Builds up a 'domination graph' by assigning each vertex to 
        its closest dominators. These dominators will be connected in the
        final graph.
    """

    h = DictGraph(nodes=domset)

    # dictionary mapping vertices from the graph to closest dominators to it
    assignment = defaultdict(set)
    # the keys of dominators should all belong to the same component, so this
    # should implicitly only operate on the vertices we care about
    for v in dominators:
        # Find the domset vertices closest to v
        sorted_doms = [dominators[v][r] for r in range(radius+1) if len(dominators[v][r])>0]
        # sorted_doms[0] contains the closest dominators (we don't care about what radius it is at)
        closest = sorted_doms[0]
        # Assign each of those closest dominating nodes to v
        for x in closest:
            assignment[v].add(x)
        # all closest dominating nodes should form a clique, since they optimally dominate a common vertex (namely, v)
        for x,y in itertools.combinations(closest,2):
            h.add_edge(x,y)

    return h, assignment

def calc_domination_graph(g, augg, domset, dominators, radius):
    """ 
        Builds up a 'domination graph' by assigning each vertex to 
        its closest dominators. These dominators will be connected in the
        final graph.
    """
    # Make the domination graph
    h, assignment = _calc_domset_graph(domset, dominators, radius)
    hcomps = h.component_index()

    #print("original assignment", assignment)

    # print("  First domgraph has {} components".format(h.num_components()))

    # We want the domination graph to be connected.  If it isn't, we need to
    # start adding edges between components that have vertices u and v
    # respectively such that uv is an edge in g
    if len(set([hcomps[v] for v in hcomps])) > 1:
        # map original graph vertices to domgraph components
        compmap = {}
        for u in dominators.keys():
            # find components of the domgraph that contain dominators of u
            candidates = [hcomps[v] for v in assignment[u]]
            #assert len(candidates) > 0
            # the components should all be the same so we can pick the first one
            assert candidates.count(candidates[0]) == len(candidates) # Make sure all the comps. agree
            compmap[u] = candidates[0]

        # look at the edges of the graph
        for u,v in g.edges():
            # domset nodes optimally dominate their neighbors
            if u in domset or v in domset:
                continue
            # vertices outside of dominators are not considered here
            if u not in dominators or v not in dominators:
                continue
            # Check whether u,v 'connects' two components of h
            # if so, add all edges between their dominators
            if compmap[u] != compmap[v]:
                for x in assignment[u]:
                    for y in assignment[v]:
                        h.add_edge(x,y)

    return h, assignment

def better_dvorak_reidl(augg,radius,comp=None):
    """ Compute a d-dominating set using Dvorak's approximation algorithm 
        for dtf-graphs (see `Structural Sparseness and Complex Networks').
        Needs a distance-d dtf augmentation augg (see rdomset() for usage). """    
    domset = set()
    infinity = float('inf')
    domdistance = defaultdict(lambda: infinity)

    # if a list of nodes in a component is supplied, we loop over that,
    # otherwise we loop over all nodes in the graph.
    if comp is None:
        nodes = augg
    else:
        nodes = comp

    #  low indegree first (reverse=False)
    vprops = [(v,augg.in_degree(v)) for v in nodes]
    vprops.sort(key=itemgetter(1),reverse=False)
    order = map(itemgetter(0),vprops)

    for v in order:
        if domdistance[v] <= radius:
            continue

        for u,r in augg.in_neighbours(v):
            domdistance[v] = min(domdistance[v],r+domdistance[u])

        if domdistance[v] <= radius:
            continue

        domset.add(v)
        domdistance[v] = 0

        for u,r in augg.in_neighbours(v):
            domdistance[u] = min(domdistance[u],r+domdistance[v])

    return domset

#
# Debugging/Test functions
#
def verify_scattered(g, sc, radius):
    """ Verifies that a vertex set sc is d-scattered in a graph g, meaning
        that every pair of vertices in sc have distance at least d. """
    sc = set(sc)
    for v in sc:
        dN = g.rneighbours(v,radius)
        if len(dN & sc) > 1: # v is in dN but nothing else should
            return False
    return True

def verify_domset(g, ds, radius):
    """ Verifies that a vertex set ds is a d-dominating set in g, meaning
        that every vertex of g lies within distance d of some vertex in ds. """
    dominated = set()
    for v in ds:
        dominated |= g.rneighbours(v,radius)
    return len(dominated) == len(g)

def verify_dominators(g, domset, dominators, radius):
    for v in g:
        dN = g.rneighbours(v,radius)
        candidates = dN & domset
        for r in range(0,radius+1):
            candidates -= dominators[v][r]
        if len(candidates) > 0:
            print(v, candidates, dominators[v])
            return False
    return True

def verify_dominators_strict(g, domset, dominators, radius):
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

def rdomset(g, radius, comp=None):
    """ Computes an d-dominating set for a given graph g
        using dtf-augmentations. """
    augg = dtf(g, radius, comp)
    domset = better_dvorak_reidl(augg, radius, comp)

    return domset, augg

def main(inputgxt, outputgxt, radius, test, domgraph):
    g, node_attrs, edge_attrs, _ = TFGraph.from_gxt(inputgxt, radius)
    g.remove_loops()

    domset, augg = rdomset(g,radius)
    print("Computed {}-domset of size {} in graph with {} vertices".format(radius,len(domset),len(g)))

    dominators = calc_dominators(augg, domset, radius)

    # Compute domset graph
    if domgraph:
        print("Computing domgraph")
        h = _calc_domset_graph(domset, dominators, radius)
        write_gxt(domgraph, h)

    # Test domset
    if test:
        print("Testing domset:", verify_domset(g,domset,radius))

    # Annotate domset
    labelds = 'ds{}'.format(d)
    for v in domset:
        node_attrs[v][labelds] = 1

    write_gxt(outputgxt, g, node_attrs, edge_attrs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='gxt input file', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('--output', help='gxt output file', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdin)
    parser.add_argument('--domgraph', help='gxt output file', nargs='?', type=argparse.FileType('w'),
                        default=None)
    parser.add_argument("radius", help="domset distance", type=int)
    parser.add_argument('--test', action='store_true')
    args = parser.parse_args()

    main(args.input, args.output, args.radius, args.test, args.domgraph)

