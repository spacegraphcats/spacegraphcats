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
        for u,__ in graph.in_neighbors(v,1):
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
