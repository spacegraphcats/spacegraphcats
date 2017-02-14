#!/usr/bin/env python3
"""Parser and writer for simple graph format."""


def _parse_line(line, split_on=' '):
    return list(map(str.strip, line.split(split_on)))


def parse(graph_file, create_vertices=None, add_edge=None):
    """Parse a graph and call provided methods with vertices and edges."""
    
    # read vertices
    num_vertices = int(graph_file.readline().strip())

    create_vertices(num_vertices);

    # read edges
    next_line = graph_file.readline()
    while len(next_line) > 1:
        parsed = _parse_line(next_line)
        add_edge(int(parsed[0]), int(parsed[1]))
        next_line = graph_file.readline()


def parse_minhash(minhash_file, add_minhash):
    """Parse minhash (.mxt) file."""
    for line in minhash_file:
        if len(line) < 2:
            continue
        parsed = _parse_line(line)
        add_minhash(parsed[0], list(map(int, map(str.strip, parsed[1:]))))

def write(graph_file, num_vertices, edges):
    """Write an edgelist into an .ext file."""
    graph_file.write('{}\n'.format(num_vertices))
    for u, v in edges:
        graph_file.write('{} {}\n'.format(u, v))
