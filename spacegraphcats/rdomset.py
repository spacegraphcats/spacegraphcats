"""Algorithms for r-dominating set computation."""
from collections import defaultdict
from spacegraphcats.graph import DictGraph


def low_degree_orientation(graph):
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


def ldo_setup(graph):
    """Setup data structure for low degree orientation."""
    n = len(graph)

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
    bins = [None for _ in graph]

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


def dtf_step(graph, dist):
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


def dtf(graph, radius):
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


def compute_domset(graph, radius):
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


def assign_to_dominators(graph, domset, radius):
    """
    Compute for each vertex the subset of domset that dominates it at
    distance radius.  Returns a double level dictionary that maps vertices in
    the original graph to integers 0 to radius to sets of vertices that
    dominate at that distance
    """
    dominated_at_radius = defaultdict(lambda: defaultdict(set))

    # Every vertex in domset is a zero-dominator of itself
    for v in domset:
        dominated_at_radius[v][0].add(v)

    # We need two passes in order to only every access the
    # in-neighbourhoods (which are guaranteed to be small).
    for _ in range(2):
        for v in graph:
            # Pull dominators from in-neighbourhood
            for u, r in graph.in_neighbors(v):
                for r2 in dominated_at_radius[u].keys():
                    domdist = r+r2
                    if domdist <= radius:
                        dominated_at_radius[v][domdist] |= \
                                        dominated_at_radius[u][r2]
            # Push dominators to in-neighbourhood
            for u, r in graph.in_neighbors(v):
                for r2 in dominated_at_radius[v].keys():
                    domdist = r+r2
                    if domdist <= radius:
                        dominated_at_radius[u][domdist] |= \
                                        dominated_at_radius[v][r2]

    # Clean up: vertices might appear at multiple distances as dominators,
    # we only want to store them for the _minimum_ distance at which they
    # act as a dominator.
    for v in dominated_at_radius:
        cumulated = set()
        for r in range(radius+1):
            dominated_at_radius[v][r] -= cumulated
            cumulated |= dominated_at_radius[v][r]

    return dominated_at_radius


def domination_graph(graph, domset, radius):
    """
    Build up a 'domination graph' by assigning each vertex to its closest
    dominators. These dominators will be connected in the final graph.
    """
    print("assigning to dominators")
    dominated_at_radius = assign_to_dominators(graph, domset, radius)
    domgraph = DictGraph(nodes=domset)

    print("computing dominating edges")
    # dictionary mapping vertices from the graph to closest dominators to it
    closest_dominators = {}
    domdistance = {v: 0 for v in dominated_at_radius}

    print("Adding edges between odd length optimal domination")
    # the keys of dominators should all belong to the same component, so this
    # should implicitly only operate on the vertices we care about
    for v in dominated_at_radius:
        # Find the domset vertices closest to v
        for r in range(radius+1):
            if len(dominated_at_radius[v][r]) > 0:
                domdistance[v] = r
                closest_dominators[v] = frozenset(dominated_at_radius[v][r])
                break
        # because they all dominate v, the closest dominators should form a
        # clique
        for x in closest_dominators[v]:
            nonneighbors = closest_dominators[v] - \
                            domgraph.in_neighbors(x, 1) - \
                            set([x])
            for y in nonneighbors:
                domgraph.add_arc(x, y)
                domgraph.add_arc(y, x)

    print("Adding edges between even length optimal domination")
    # Two vertices in the dominating graph should be adjacent if they
    # optimally dominate a common vertex.  Note that this only includes
    # dominators with odd length paths between them because even length paths
    # have no midpoint that would be dominated at the same distance from both
    # directions.  So we need to expand our criteria to include dominators
    # that optimally dominate neighbors at the same distance.  Since the
    # endpoints of each edge have to be optimally dominated at the same radius
    # or at radii that differ by 1, it is sufficient to loop through all edges
    # of the original graph and turn the set of dominators into a clique.
    # This will also catch and connect the dominators that were themselves
    # adjacent in the original graph.
    for v in graph:
        for u in graph.in_neighbors(v, 1):
            if domdistance[u] != domdistance[v]:
                continue
            # since the edges within closest_dominators have been added above,
            # we only need to look at the symmetric difference
            difference = closest_dominators[u] ^ closest_dominators[v]
            for x in difference:
                nonneighbors = difference - \
                                domgraph.in_neighbors(x, 1) - \
                                set([x])
                for y in nonneighbors:
                    domgraph.add_arc(x, y)
                    domgraph.add_arc(y, x)

    return domgraph, closest_dominators


def rdomset(graph, radius):
    """Compute a dominating set of graph at radius radius."""
    dtf(graph, radius)
    domset = compute_domset(graph, radius)
    return domset
