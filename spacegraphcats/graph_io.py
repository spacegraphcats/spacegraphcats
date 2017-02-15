from spacegraphcats.graph_parser import parse, write
from spacegraphcats.graph import Graph

def read_from_gxt(gxtfile, radius: int):
    """
    Read an gxt file and returns a graph with the edges in the mxt file.
    """
    g = None  # type: Graph

    def create_graph(num_nodes):
        nonlocal g
        g = Graph(num_nodes, radius)

    def add_edge(u, v):
        g.add_arc(u, v)

    parse(gxtfile, create_graph, add_edge)
    
    return g

def write_to_gxt(gxtfile, graph, weight: int = None):
    """
    Write graph to gxt file.  If a weight is provided, only the arcs with a particular weight are stored.
    """
    write(gxtfile, len(graph), graph.arcs(weight))
