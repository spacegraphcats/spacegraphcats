from spacegraphcats.graph_parser import parse, write
from spacegraphcats.graph import Graph

def read_from_gxt(gxtfile, radius: int, directed: bool):
    """
    Read an gxt file and returns a graph with each edge in the mxt file.
    When directed is set to False, there will be two edges for each edge in the file.
    """
    graph = None  # type: Graph

    def create_graph(num_nodes):
        nonlocal graph
        graph = Graph(num_nodes, radius)

    def add_edge(u, v):
        graph.add_arc(u, v)
        if not directed:
            graph.add_arc(v, u)

    parse(gxtfile, create_graph, add_edge)
    
    return graph

def write_to_gxt(gxtfile, graph, weight: int = None):
    """
    Write graph to gxt file as an undorected graph.
    If a weight is provided,only the arcs with a particular weight are stored.
    """
    write(gxtfile, len(graph), [(x[0], x[1]) for x in graph.arcs(weight)])
