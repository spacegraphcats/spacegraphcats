from collections import defaultdict

from graph import DictGraph

def low_degree_orientation(graph, comp=None):
    """ 
    Computes a low-in-degree orientation of a graph component in place
    by iteratively removing a vertex of mimimum degree and orienting
    the edges towards it. 
    Precondition:  every edge in the component has a corresponding arc in the
    anti-parallel direction (i.e. uv and vu are in the graph)
    """

    # if a list of nodes in a component is supplied, we loop over that,
    # otherwise we loop over all nodes in the graph.
    if comp is None:
        nodes = range(len(graph))
    else:
        nodes = comp

    # precompute the degrees of each vertex and make a bidirectional lookup
    degdict = {}
    buckets = defaultdict(set)
    for v in nodes:
        d = graph.in_degree(v)
        degdict[v] = d
        buckets[d].add(v)

    # keep track of "removed" vertices
    removed = set()
    # run the loop once per vertex
    for _ in nodes:
        # find the minimum degree of all vertices
        d = 0
        while len(buckets[d]) == 0:
            d += 1
        # grab a vertex of minimum degree
        v = buckets[d].pop()

        # decrement the degrees of the in neighbors not yet removed and orient
        # edges towards in neighbors already removed
        for u in graph.in_neighbors(v,1):
            # if we've removed u, we orient the arc towards u by deleting uv
            if u in removed:
                graph.remove_arc(u,v)
            # otherwise, the effective degree of u should drop by 1
            else:
                d = degdict[u]
                buckets[d].remove(u)
                buckets[d-1].add(u)
                degdict[u] -= 1

        removed.add(v)

def dtf_step(graph, dist, comp=None):
    """ 
    Computes the d-th dtf-augmentation from a graph 
    (where d is provided by the argument dist). The input graph
    must be a (d-1)-th dtf-augmentation. See dtf() for usage. 
    This function adds arcs to graph.
    """

    fratGraph = DictGraph() # Records fraternal edges, must be oriented at the end
    newTrans = set() # Records transitive arcs

    # if a list of nodes in a component is supplied, we loop over that,
    # otherwise we loop over all nodes in the graph.
    if comp is None:
        nodes = range(len(graph))
    else:
        nodes = comp

    # pick out the transitive and fraternal pairs from v
    for v in nodes:
        for x, y in graph.transitive_pairs(v, dist):
            #assert x != y
            newTrans[(x, y)] = dist
        for x, y in graph.fraternal_pairs(v, dist):
            #assert x != y
            fratGraph.add_arc(x, y)
            fratGraph.add_arc(y, x)

    # Add transitive arcs to graph
    for (s, t) in newTrans:
        #assert s != t
        augg.add_arc(s, t, dist)
        # if s,t is both a transitive and fraternal pair, don't add it twice!
        fratGraph.remove_arc(s,t)
        fratGraph.remove_arc(t,s)

    # Orient fraternal edges and add them to the grah
    low_degree_orientation(fratGraph)

    for s, t in fratDigraph.arcs(1):
        #assert s != t
        graph.add_arc(s,t,dist)

def dtf(graph, radius, comp=None):
    """
    Computes dft-augmentations of a graph.  The number of augmentations is
    given by the argument radius.
    Postcondition:
        If the distance between each pair of vertices u,v is at most radius in
        the original graph, uv or vu will be an arc with weight equal to that
        distance.
    """
    # the 1st "augmentation" is simply acyclically orienting the edges
    low_degree_orientation(graph ,comp=comp)

    # keep track of whether we are adding edges so we can quit early
    num_arcs = graph.num_arcs()
    changed = True
    d = 2
    while changed and d <= radius:
        # shortcut paths of length d
        dtf_step(graph, d, comp)

        # Small optimization: if no new arcs have been added we can stop.
        curr_arcs = graph.num_arcs() # This costs a bit so we store it
        changed = num_arcs < curr_arcs
        num_arcs = curr_arcs
        d += 1

def compute_domset(graph,radius,comp=None):
    """ Compute a d-dominating set using Dvorak's approximation algorithm 
        for dtf-graphs (see `Structural Sparseness and Complex Networks').
        graph needs a distance-d dtf augmentation (see rdomset() for usage). """    
    domset = set()
    infinity = float('inf')
    # minimum distance to a dominating vertex, obviously infinite at start
    domdistance = defaultdict(lambda: infinity)

    # if a list of nodes in a component is supplied, we loop over that,
    # otherwise we loop over all nodes in the graph.
    if comp is None:
        nodes = graph
    else:
        nodes = comp

    # Sort the vertices by indegree so we take fewer vertices
    order = sorted([v for v in nodes],key=lambda x:graph.in_degree(x))
    # vprops = [(v,graph.in_degree(v)) for v in nodes]
    # vprops.sort(key=itemgetter(1),reverse=False)
    # order = map(itemgetter(0),vprops)

    for v in order:
        # if v is already dominated at radius, no need to work
        if domdistance[v] <= radius:
            continue

        # look at the in neighbors to update the distance
        for u,r in graph.in_neighbors(v):
            domdistance[v] = min(domdistance[v],r+domdistance[u])

        # if v is dominated at radius now, keep going
        if domdistance[v] <= radius:
            continue
        # otherwise put v in the dominating set
        domset.add(v)
        domdistance[v] = 0

        # update distances of neighbors of v if v is closer
        for u,r in graph.in_neighbors(v):
            domdistance[u] = min(domdistance[u],r)

    return domset

def assign_to_dominators(graph, domset, radius,comp=None):
    """
    Computes for each vertex the subset of domset that dominates it at distance radius
    Returns a double level dictionary that maps vertices in the original graph to
    integers 0 to radius to sets of vertices that dominate at that distance
    """
    dominated_at_radius = defaultdict(lambda: defaultdict(set))

    if comp is None:
        nodes = graph
    else:
        nodes = comp

    # Every vertex in domset is a zero-dominator of itself
    for v in domset:
        dominated_at_radius[v][0].add(v)

    # We need two passes in order to only every access the
    # in-neighbourhoods (which are guaranteed to be small).
    for _ in range(2):
        for v in nodes:
            # Pull dominators from in-neighbourhood
            for u, r in graph.in_neighbors(v):
                for r2 in dominated_at_radius[u].keys():
                    domdist = r+r2
                    if domdist <= radius:
                        dominated_at_radius[v][domdist] |= dominated_at_radius[u][r2]
            # Push dominators to in-neighbourhood
            for u, r in graph.in_neighbors(v):
                for r2 in dominated_at_radius[v].keys():
                    domdist = r+r2
                    if domdist <= radius:
                        dominated_at_radius[u][domdist] |= dominated_at_radius[v][r2]

    # Clean up: vertices might appear at multiple distances as dominators,
    # we only want to store them for the _minimum_ distance at which they
    # act as a dominator.
    for v in dominated_at_radius:
        cumulated = set()
        for r in range(radius+1):
            dominated_at_radius[v][r] -= cumulated
            cumulated |= dominated_at_radius[v][r]

    return dominated_at_radius

def domination_graph(graph, domset, dominated_at_radius):
    """ 
        Builds up a 'domination graph' by assigning each vertex to 
        its closest dominators. These dominators will be connected in the
        final graph.
        Precondition:
            The keys of dominated_at_radius are exactly one connected 
            component of graph
    """

    domgraph = DictGraph(nodes=domset)

    # dictionary mapping vertices from the graph to closest dominators to it
    closest_dominators = {v:list() for v in dominated_at_radius}
    # the keys of dominators should all belong to the same component, so this
    # should implicitly only operate on the vertices we care about
    for v in dominated_at_radius:
        # Find the domset vertices closest to v
        sorted_doms = [dominated_at_radius[v][r] for r in range(radius+1) if len(dominated_at_radius[v][r])>0]
        # sorted_doms[0] contains the closest dominators (we don't care about 
        # what radius it is at)
        closest = sorted_doms[0]
        # Assign each of those closest dominating nodes to v
        for x in closest:
            closest_dominators[v].append(x)
        # all closest dominating nodes should form a clique, since they 
        # optimally dominate a common vertex (namely, v)
        for x,y in itertools.combinations(closest,2):
            domgraph.add_arc(x,y)
            domgraph.add_arc(y,x)
    # all adjacent vertices in the dominating set should also be adjacent in 
    # the dominating graph.  If two dominators are adjacent, they probably 
    # optimally a common vertex and are already connected anyway, but we want 
    # to make sure we get all of them
    for v in domset:
        adj_doms = dominated_at_radius[v][1]
        for u in adj_doms:
            domgraph.add_arc(u,v)
            domgraph.add_arc(v,u)

    # ensure domgraph is connected
    make_connected(domgraph, domset, closest_dominators, graph)

    return domgraph, closest_dominators

def make_connected(domgraph, domset, closest_dominators, graph):
    """
    Makes a domination graph connected.  If it isn't already connected, we 
    need to start adding edges between components that have vertices u and v 
    respectively such that uv is an edge in graph.  Note that this will only 
    occur if u and v are optimally dominated by x and y respectively at the 
    same distance.
    """
    # map vertices to indices of components
    dom_components = domgraph.component_index()
    num_comps = len(set([dom_components[v] for v in domset]))
    # don't do anything it it's already connected!
    if num_comps == 1:
        return
    # find the components the non-dominating vertices' dominators belong to
    nondom_components = {}
    for u in closest_dominators.keys():
        # since all of u's dominators are in the same component, we can pick 
        # an arbitrary one to look up
        nondom_components[u] = dom_components[closest_dominators[u][0]]

    # look at arcs in graph and determine whether they bridge components
    for u in closest_dominators.keys():
        # the domset vertices optimally dominate their neighbors, so they 
        # can't help us merge components
        if u in domset:
            continue
        for v in graph.in_neighbors(u,1):
            # if this edge bridges components, connect all the optimal 
            # dominators
            if nondom_components[u] != nondom_components[v]:
                for x in closest_dominators[u]:
                    for y in closest_dominators[v]:
                        domgraph.add_arc(x,y)
                        domgraph.add_arc(y,x)

def rdomset(graph, radius, comp=None):
    dtf(graph, radius, comp)
    domset = compute_domset(graph, radius, comp)


