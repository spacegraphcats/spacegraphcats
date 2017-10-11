"""Reading and writing functions."""
from .graph_parser import parse, write
from .graph import Graph, DictGraph


def read_from_gxt(gxtfile, radius: int, directed: bool, sequential=True):
    """
    Read an gxt file and returns a graph with each edge in the mxt file.

    When directed is set to False, there will be two edges for each edge in
    the file.
    """
    graph = None  # type: Graph

    def create_graph(num_nodes):
        nonlocal graph
        if sequential:
            print("Sequential graph with {} nodes ".format(num_nodes))
            graph = Graph(num_nodes, radius)
        else:
            graph = DictGraph(r=radius)

    def add_edge(u, v):
        if u != v:
            if not sequential:
                if u not in graph:
                    graph.add_node(u)
                if v not in graph:
                    graph.add_node(v)
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
    write(gxtfile, len(graph), map(lambda x: (x[0], x[1]), graph.arcs(weight)))
