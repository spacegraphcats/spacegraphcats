"""Algorithms for r-dominating set computation."""
from collections import defaultdict
from .graph import Graph, DictGraph

from typing import List, Set, Dict, Union

from sortedcontainers import SortedSet, SortedDict


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
    (
        bins,
        # pointers to the first vertex with a given degree
        bin_starts,
        # remaining degree of the vertex
        degrees,
        # pointer to the location of a vertex in bins
        location,
    ) = ldo_setup(graph)

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
            print("removed {} of {} nodes\r".format(curr + 1, n), end="")
            checkpoint += n // 100
    print("removed {} of {} nodes".format(curr + 1, n))


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
    degree_counts = [0 for i in range(max_deg + 1)]
    for v in graph:
        d = degrees[v]
        degree_counts[d] += 1
    # assign the cutoffs of bins
    bin_starts = [sum(degree_counts[:i]) for i in range(max_deg + 1)]
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
            print("bucketed {} of {} nodes\r".format(i + 1, n), end="")
            checkpoint += n // 100
    del bin_ptrs
    print("bucketed {} of {} nodes".format(i + 1, n))
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
    print("added {} fraternal edges".format(fratGraph.num_arcs() // 2))

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
    domset = SortedSet()
    infinity = float("inf")
    # minimum distance to a dominating vertex, obviously infinite at start
    domdistance = defaultdict(lambda: infinity)  # type: Dict[int, float]
    # counter that keeps track of how many neighbors have made it into the
    # domset
    domcounter = defaultdict(int)  # type: Dict[int, int]
    # cutoff for how many times a vertex needs to have its neighbors added to
    # the domset before it does.  We choose radius^2 as a convenient "large"
    # number
    c = (2 * radius) ** 2

    # Sort the vertices by indegree so we take fewer vertices
    order = sorted([v for v in graph], key=lambda x: graph.in_degree(x), reverse=True)
    # vprops = [(v,graph.in_degree(v)) for v in nodes]
    # vprops.sort(key=itemgetter(1),reverse=False)
    # order = map(itemgetter(0),vprops)

    for v in order:
        # look at the in neighbors to update the distance
        for r in range(1, radius + 1):
            for u in graph.in_neighbors(v, r):
                domdistance[v] = min(domdistance[v], r + domdistance[u])

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


def domination_graph(graph: Graph, domset: Set[int], radius: int):
    """
    Build up a 'domination graph' by assigning each vertex to one of its
    closest dominators, then preserving edges between vertices assigned to
    different dominators.  This means the domination graph is a shallow minor
    (at radius radius) of the original graph.
    """
    print("assigning to dominators")
    domgraph = DictGraph(nodes=SortedSet(domset))

    # We need to assign each vertex to a unique closest dominator.  A value of
    # -1 indicates this vertex has not yet been assigned to a dominator.
    assigned_dominator = {v: -1 for v in graph}  # type: Dict[int, int]
    # domset vertices are closest to themselves
    for v in domset:
        assigned_dominator[v] = v

    # If we compute all closest dominators and arbitrarily pick a unique one
    # from the list, the vertices assigned to a particular dominator
    # (including the dominator itself) might not induce a connected subgraph.
    # Thus we need to propagate the domination assignments outwards so a
    # vertex at distance d from a dominator has to wait until its neighbors at
    # distance d-1 are assigned in order to pick an appropriate assignment (a
    # dominator already found in its neighborhood).
    for _ in range(radius):
        for v in graph:
            label = assigned_dominator[v]
            if label >= 0:
                # Push value to unassigned in-neighbors
                for u in graph.in_neighbors(v, 1):
                    if assigned_dominator[u] < 0:
                        assigned_dominator[u] = label
            else:
                # Pull value
                for u in graph.in_neighbors(v, 1):
                    # TODO: could instead use some min/max voting here
                    if assigned_dominator[u] >= 0:
                        assigned_dominator[v] = assigned_dominator[u]
                        break
    # Sanity-check
    for v, value in enumerate(assigned_dominator):
        assert value >= 0

    # Compute domgraph edges
    print("computing domgraph edges")
    for u, v in graph.arcs(1):
        du, dv = assigned_dominator[u], assigned_dominator[v]
        assert du in domset and dv in domset
        if du != dv:
            domgraph.add_arc(du, dv)
            domgraph.add_arc(dv, du)

    # map of dominators to vertices they dominate
    dominated = SortedDict(
        {x: SortedSet() for x in domset}
    )  # type: Dict[int, Set[int]]
    for v, x in assigned_dominator.items():
        dominated[x].add(v)

    domgraph.remove_isolates()

    return domgraph, dominated


def rdomset(graph: Graph, radius: int):
    """Compute a dominating set of graph at radius radius."""
    dtf(graph, radius)
    domset = compute_domset(graph, radius)
    return domset
