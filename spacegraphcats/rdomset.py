"""Algorithms for r-dominating set computation."""
from collections import defaultdict
from .graph import Graph, DictGraph

from typing import List, Set, Dict, Any, Union


def low_degree_orientation(graph: Graph):
    """
    Compute a low-in-degree orientation of a graph.

    This occurs in place by iteratively removing a vertex of mimimum degree
    and orienting the edges towards it.

    Precondition:  every edge in the component has a corresponding arc in the
    anti-parallel direction (i.e. uv and vu are in the graph)
    """
    # number of vertices (needed for print statements)
    n = len(graph)
    # most elegant way to handle possibly empty graph (from the frat graph)
    if n == 0:
        return
    """
    compute necessary data structures for low degree orientation
    """
    # array binning the vertices by remaining degree
    (bins,
        # pointers to the first vertex with a given degree
        bin_starts,
        # remaining degree of the vertex
        degrees,
        # pointer to the location of a vertex in bins
        location) = ldo_setup(graph)

    checkpoint = 0
    # run the loop once per vertex
    for curr in range(n):
        # curr points the vertex of minimum degree
        v = bins[curr]
        d_v = degrees[v]
        # "move" v into bin 0 if it isn't there
        if d_v > 0:
            for i in range(d_v, 0, -1):
                bin_starts[i] += 1
        degrees[v] = 0
        # decrement the degrees of the in neighbors not yet removed and orient
        # edges towards in neighbors already removed
        inbrs = list(graph.in_neighbors(v, 1))
        for u in inbrs:
            d_u = degrees[u]
            loc_u = location[u]
            # if we've removed u, we orient the arc towards u by deleting uv
            if location[u] < location[v]:
                graph.remove_arc(u, v)
            # otherwise, the effective degree of u should drop by 1
            else:
                # swap u with w, the first vertex with the same degree
                # find where w is
                loc_w = bin_starts[d_u]
                w = bins[loc_w]
                # swap their positions
                if w != u:
                    bins[loc_u] = w
                    bins[loc_w] = u
                    location[w] = loc_u
                    location[u] = loc_w
                # move the bin start one place over
                bin_starts[d_u] += 1
                # decrement u's degree
                degrees[u] = d_u - 1
        if curr == checkpoint:
            print("removed {} of {} nodes\r".format(curr+1, n), end="")
            checkpoint += n//100
    print("removed {} of {} nodes".format(curr+1, n))


def ldo_setup(graph: Graph):
    """Setup data structure for low degree orientation."""
    n = len(graph)

    degrees = None  # type: Union[Dict[int, int], List[int]]
    location = None  # type: Union[Dict[int, int], List[int]]
    max_deg = 0

    # hack-y way to know whether our location and degree lookups should be
    # lists or dictionaries
    if isinstance(graph, DictGraph):
        # degree lookup
        degrees = {v: graph.in_degree(v) for v in graph}
        # pointer to place in vertex ordering
        location = {v: None for v in graph}
        max_deg = max(degrees.values())

    else:
        degrees = [graph.in_degree(v) for v in graph]
        location = [None for _ in graph]
        max_deg = max(degrees)

    # precompute the degrees of each vertex and make a bidirectional lookup
    degree_counts = [0 for i in range(max_deg+1)]
    for v in graph:
        d = degrees[v]
        degree_counts[d] += 1
    # assign the cutoffs of bins
    bin_starts = [sum(degree_counts[:i]) for i in range(max_deg+1)]
    del degree_counts
    bin_ptrs = list(bin_starts)
    bins = [None for _ in graph]  # type: List[int]

    # assign the vertices to bins
    checkpoint = 0
    for i, v in enumerate(graph):
        loc = bin_ptrs[degrees[v]]
        bins[loc] = v
        location[v] = loc
        bin_ptrs[degrees[v]] += 1
        if v == checkpoint:
            print("bucketed {} of {} nodes\r".format(i+1, n), end="")
            checkpoint += n//100
    del bin_ptrs
    print("bucketed {} of {} nodes".format(i+1, n))
    return bins, bin_starts, degrees, location


def dtf_step(graph: Graph, dist):
    """
    Compute the d-th dtf-augmentation from a graph.

    d is provided by the argument dist. The input graph
    must be a (d-1)-th dtf-augmentation. See dtf() for usage.
    This function adds arcs to graph.
    """
    fratGraph = DictGraph()
    # Records fraternal edges, must be oriented at the end

    trans_pairs = 0
    # pick out the transitive pairs from v and add them as new edges
    for v in graph:
        for x, y in graph.transitive_pairs(v, dist):
            graph.add_arc(x, y, dist)
            trans_pairs += 1
    print("added {} transitive edges".format(trans_pairs))
    # pick out the fraternal pairs from v and store them.  We do this after
    # adding transitive edges to guarantee that no fraternal edge conflicts
    # with a transitive edge
    for v in graph:
        for x, y in graph.fraternal_pairs(v, dist):
            # assert x != y
            fratGraph.add_node(x)
            fratGraph.add_node(y)
            fratGraph.add_arc(x, y)
            fratGraph.add_arc(y, x)
    print("added {} fraternal edges".format(fratGraph.num_arcs()//2))

    # Orient fraternal edges and add them to the graph
    low_degree_orientation(fratGraph)

    for s, t in fratGraph.arcs(1):
        # assert s != t
        graph.add_arc(s, t, dist)


def dtf(graph: Graph, radius: int):
    """
    Compute dft-augmentations of a graph.

    The number of augmentations is given by the argument radius.
    Postcondition:
        If the distance between each pair of vertices u,v is at most radius in
        the original graph, uv or vu will be an arc with weight equal to that
        distance.
    """
    # the 1st "augmentation" is simply acyclically orienting the edges
    print("Computing low degree orientation (step 1)")
    low_degree_orientation(graph)

    # keep track of whether we are adding edges so we can quit early
    # num_arcs = graph.num_arcs()
    changed = True
    d = 2
    while changed and d <= radius:
        print("Computing step {}".format(d))
        # shortcut paths of length d
        dtf_step(graph, d)

        # Small optimization: if no new arcs have been added we can stop.
        # curr_arcs = graph.num_arcs() # This costs a bit so we store it
        # changed = num_arcs < curr_arcs
        # num_arcs = curr_arcs
        d += 1


def compute_domset(graph: Graph, radius: int):
    """
    Compute a d-dominating set using Dvorak's approximation algorithm
    for dtf-graphs (see `Structural Sparseness and Complex Networks').
    Graph needs a distance-d dtf augmentation (see rdomset() for usage).
    """
    domset = set()
    infinity = float('inf')
    # minimum distance to a dominating vertex, obviously infinite at start
    domdistance = defaultdict(lambda: infinity)
    # counter that keeps track of how many neighbors have made it into the
    # domset
    domcounter = defaultdict(int)
    # cutoff for how many times a vertex needs to have its neighbors added to
    # the domset before it does.  We choose radius^2 as a convenient "large"
    # number
    c = (2*radius)**2

    # Sort the vertices by indegree so we take fewer vertices
    order = sorted([v for v in graph], key=lambda x: graph.in_degree(x),
                   reverse=True)
    # vprops = [(v,graph.in_degree(v)) for v in nodes]
    # vprops.sort(key=itemgetter(1),reverse=False)
    # order = map(itemgetter(0),vprops)

    for v in order:
        # look at the in neighbors to update the distance
        for r in range(1, radius + 1):
            for u in graph.in_neighbors(v, r):
                domdistance[v] = min(domdistance[v], r+domdistance[u])

        # if v is already dominated at radius, no need to work
        if domdistance[v] <= radius:
            continue

        # if v is not dominated at radius, put v in the dominating set
        domset.add(v)
        domdistance[v] = 0

        # update distances of neighbors of v if v is closer if u has had too
        # many of its neighbors taken into the domset, include it too.
        for r in range(1, graph.radius + 1):
            for u in graph.in_neighbors(v, r):
                domcounter[u] += 1
                domdistance[u] = min(domdistance[u], r)
                if domcounter[u] > c and u not in domset:
                    # add u to domset
                    domset.add(u)
                    domdistance[u] = 0
                    for x, rx in graph.in_neighbors(u):
                        domdistance[x] = min(domdistance[x], rx)
                # only need to update domdistance if u didn't get added

    return domset


def assign_to_dominators(graph: Graph, domset: Set[int], radius: int):
    """
    Compute the closest dominators to each vertex.
    
    Return values:
        closest_dominators:  dictionary mapping vertices in the original graph
                             to their closest_dominators
        domdistance:  dictionary mapping vertices to the distance to their
                      closest dominators
    """
    closest_dominators = {v:set() for v in graph} # type: Dict[int, Set[int]]
    domdistance = {v:radius+1 for v in graph} # type: Dict[int, int]

    # Every vertex in domset is a zero-dominator of itself
    for v in domset:
        closest_dominators[v] = set([v])
        domdistance[v] = 0

    # We need two passes in order to only every access the
    # in-neighbourhoods (which are guaranteed to be small).
    for _ in range(2):
        for v in graph:
            # Pull dominators from in-neighbourhood
            for u, r in graph.in_neighbors(v):
                dist = r + domdistance[u]
                # we need to reset v's dominators when the path through u
                # yields a closer dominator
                if dist < domdistance[v]:
                    closest_dominators[v] = set()
                    domdistance[v] = dist
                # add u's dominators to v's
                if dist == domdistance[v]:
                    closest_dominators[v] |= closest_dominators[u]
            # Push dominators to in-neighbourhood
            for u, r in graph.in_neighbors(v):
                dist = r + domdistance[v]
                # we need to reset u's dominators when the path through v
                # yields a closer dominator
                if dist < domdistance[u]:
                    closest_dominators[u] = set()
                    domdistance[u] = dist
                # add v's dominators to u's
                if dist == domdistance[u]:
                    closest_dominators[u] |= closest_dominators[v]

    return closest_dominators, domdistance


def domination_graph(graph: Graph, domset: Set[int], radius: int):
    """
    Build up a 'domination graph' by assigning each vertex to its closest
    dominators. These dominators will be connected in the final graph.
    """
    print("assigning to dominators")
    closest_dominators, domdistance = assign_to_dominators(graph, domset, radius)
    domgraph = DictGraph(nodes=domset)

    print("computing dominating edges")
    # map of dominators to vertices they dominate
    dominated = {x:set() for x in domset} # type: Dict[int, Set[int]]
    for v, doms in closest_dominators.items():
        for x in doms:
            dominated[x].add(v)

    # Add edges to the domgraph.  Two dominators (x,y) should be adjacent if
    # they optimally dominate a common vertex or if there are adjacent vertices
    # (u,v) optimally dominated by (x,y), respectively, at the same radius.
    # If v is optimally dominated by x at radius r, all of its neighbors are
    # dominated at radius r-1, r, r+1.  If u is a neighbor of v not optimally 
    # dominated by x, then either:
    #    1) u is dominated at distance r OR 
    #    2) u is dominated at distance r-1 and v's optimal dominators are a
    #       superset of u's
    # In either case x should be adjacent to all dominators of u.
    for x in domset:
        # neighbors of vertices dominated by x that are not dominated by x
        domination_boundary = set() # type: Set[int]
        # vertices that need to be made adjacent to x
        new_dom_neighbors = set() # type: Set[int]
        # get all the boundary vertices.  We remove the vertices dominated by x
        # at the end to only do it once
        for v in dominated[x]:
            domination_boundary |= graph.in_neighbors(v, 1)
        domination_boundary -= dominated[x]
        # get the dominators of the boundary.
        for v in domination_boundary:
            new_dom_neighbors |= closest_dominators[v]
        # take out existing neighbors
        new_dom_neighbors -= domgraph.in_neighbors(x, 1)
        # add edges
        for y in new_dom_neighbors:
            domgraph.add_arc(x, y)
            domgraph.add_arc(y, x)
    
    return domgraph, closest_dominators


def rdomset(graph: Graph, radius: int):
    """Compute a dominating set of graph at radius radius."""
    dtf(graph, radius)
    domset = compute_domset(graph, radius)
    return domset
