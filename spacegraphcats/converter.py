"""Convert graph with nonconsecutive ids to have consecutive ids."""
import sys

from old.spacegraphcats.graph_parser import parse
from spacegraphcats.graph_parser import write


def main():
    """Apply to graph."""
    edges = []
    num_vertices = 0

    def add_edge(u, v, *args):
        edges.append((u, v))

    def add_vertex(u, *args):
        nonlocal num_vertices
        num_vertices += 1

    parse(sys.stdin, add_vertex, add_edge, consecutive_ids=True)
    write(sys.stdout, num_vertices, edges)


if __name__ == '__main__':
    main()
