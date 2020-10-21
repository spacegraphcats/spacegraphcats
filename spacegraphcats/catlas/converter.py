#!/usr/bin/env python3

"""Convert graph with nonconsecutive ids to have consecutive ids."""
import sys

from .graph_parser import write


class IdentityHash:
    def __init__(self):
        pass

    def __getitem__(self, item):
        return item


def _parse_line(line):
    return list(map(str.strip, line.split(",")))


def parse(graph_file, add_vertex=None, add_edge=None, consecutive_ids=False):
    """Parser for (old) simple graph format.

    Parse a graph and call provided methods with vertices and edges."""
    # read vertices
    vertex_attributes = _parse_line(graph_file.readline())[2:]

    # consecutive id to original id
    if consecutive_ids:
        id_map = []
    else:
        id_map = IdentityHash()
    # original id to consecutive id
    id_map_reverse = {}

    def _get_consecutive_id(id):
        if not consecutive_ids:
            return id

        if id in id_map_reverse:
            return id_map_reverse[id]
        else:
            consecutive_id = len(id_map)
            id_map_reverse[id] = consecutive_id
            id_map.append(id)
            return consecutive_id

    next_line = graph_file.readline()
    while len(next_line) > 1:
        if add_vertex is not None:
            parsed = _parse_line(next_line)
            add_vertex(
                _get_consecutive_id(int(parsed[0])),
                int(parsed[1]),
                vertex_attributes,
                parsed[2:],
            )
        next_line = graph_file.readline()

    if add_edge is None:
        # we won't be doing anything with the edges anyway
        return id_map

    # read edges
    edge_attributes = _parse_line(graph_file.readline())[2:]

    next_line = graph_file.readline()
    while len(next_line) > 1:
        parsed = _parse_line(next_line)
        add_edge(
            _get_consecutive_id(int(parsed[0])),
            _get_consecutive_id(int(parsed[1])),
            edge_attributes,
            parsed[2:],
        )
        next_line = graph_file.readline()

    return id_map


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


if __name__ == "__main__":
    main()
