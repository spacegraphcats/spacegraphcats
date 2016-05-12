#!/usr/bin/env python3
"""Parser for simple graph format."""


def _parse_line(line):
    return list(map(str.strip, line.split(',')))


def _mapstr(l):
    return list(map(str, l))


def parse(file, add_vertex, add_edge):
    """Parse a graph and call provided methods with vertices and edges."""
    # read vertices
    vertex_attributes = _parse_line(file.readline())[2:]

    next_line = file.readline()
    while len(next_line) > 1:
        parsed = _parse_line(next_line)
        add_vertex(int(parsed[0]), int(parsed[1]),
                   vertex_attributes, parsed[2:])
        next_line = file.readline()

    # read edges
    edge_attributes = _parse_line(file.readline())[2:]

    next_line = file.readline()
    while len(next_line) > 1:
        parsed = _parse_line(next_line)
        add_edge(int(parsed[0]), int(parsed[1]),
                 edge_attributes, parsed[2:])
        next_line = file.readline()


def parse_minhash(file, add_minhash):
    """Parse minhash (.mxt) file."""
    for line in file:
        if len(line) < 2:
            continue
        parsed = _parse_line(line)
        add_minhash(parsed[0], list(map(int,map(str.strip, parsed[1].split()))))

def _parse_edgelist(file, add_edge):
    """Parse and edgelist (.ext) file."""
    for line in file:
        parsed = _parse_line(line)
        add_edge(parsed[0], parsed[1])

def write_edgelist(file, edges):
    """Write an edgelist into an .ext file."""
    for u,v in edges:
        file.write('{},{}\n'.format(u,v))

class Writer:
    """Writer for the gxt graph format.

    You need to either pass the vertex and edge attributes (the names)
    when you initialize the writer or when you add a vertex or edge.
    """

    def __init__(self, file, vertex_attributes=None, edge_attributes=None):
        """Initialize graph writer."""
        self.file = file
        self.vertex_header_written = False
        self.edge_header_written = False

        if vertex_attributes is not None:
            self.vertex_header = ','.join(['id', 'size'] + vertex_attributes)

        if edge_attributes is not None:
            self.edge_header = ','.join(['src', 'dest'] + edge_attributes)

    def add_vertex(self, id, size, attribute_values=[], vertex_attributes=None):
        """Add a vertex to the output. Don't add edges after adding nodes."""
        assert not self.edge_header_written

        if not hasattr(self, 'vertex_header'):
            self.vertex_header = ','.join(['id', 'size'] + vertex_attributes)

        if not self.vertex_header_written:
            self.file.write(self.vertex_header + '\n')
            self.vertex_header_written = True
        self.file.write(','.join(_mapstr([id, size] + attribute_values)) + '\n')

    def add_edge(self, src, dest, attribute_values=[], edge_attributes=None):
        """Add an edge to the output. Add all the nodes before you add edges."""
        assert self.vertex_header_written

        if not hasattr(self, 'edge_header'):
            self.edge_header = ','.join(['src', 'dest'] + edge_attributes)

        if not self.edge_header_written:
            self.file.write('\n')
            self.file.write(self.edge_header + '\n')
            self.edge_header_written = True
        self.file.write(','.join(_mapstr([src, dest] + attribute_values)) + '\n')

    def done(self):
        """Call when done."""
        pass


class GmlWriter:
    """Similar to the writer for gxt above but for gml."""

    def __init__(self, file, vertex_attributes=None, edge_attributes=None, directed=False):
        """Initialize graph writer."""
        self.file = file

        if vertex_attributes is not None:
            self.vertex_attributes = vertex_attributes
        if edge_attributes is not None:
            self.edge_attributes = edge_attributes

        if directed:
            self._write('graph [\n   directed 1\n')
        else:
            self._write('graph [\n   directed 0\n')

    def _write(self, string):
        self.file.write(string)

    def _quote(self, value):
        if isinstance(value, str):
            return '"{}"'.format(value)
        return value

    def add_vertex(self, id, size, attribute_values=[], vertex_attributes=None):
        """Add a vertex to the output."""
        if not hasattr(self, 'vertex_attributes'):
            self.vertex_attributes = vertex_attributes

        self._write('  node [\n')
        self._write('    id {}\n'.format(id))
        self._write('    size {}\n'.format(size))
        for k, v in zip(self.vertex_attributes, attribute_values):
            self._write('    {} {}\n'.format(k, self._quote(v)))
        self._write('  ]\n')

    def add_edge(self, src, dest, attribute_values=[], edge_attributes=None):
        """Add an edge to the output."""
        if not hasattr(self, 'edge_attributes'):
            self.edge_attributes = edge_attributes

        self._write('  edge [\n')
        self._write('    source {}\n'.format(src))
        self._write('    target {}\n'.format(dest))
        for k, v in zip(self.edge_attributes, attribute_values):
            self._write('    {} {}\n'.format(k, self._quote(v)))
        self._write('  ]\n')

    def done(self):
        """Call when done."""
        self._write(']\n')


class DotWriter:
    """Similar to the writer for gxt above but for dot."""

    def __init__(self, file):
        """Initialize graph writer."""
        self.file = file

        self._write('graph G {\n')

    def _write(self, string):
        self.file.write(string)

    def add_vertex(self, id):
        """Add a vertex to the output."""
        self._write('  {};\n'.format(id))

    def add_edge(self, src, dest):
        """Add an edge to the output."""
        self._write('  {} -- {};\n'.format(src, dest))

    def done(self):
        """Call when done."""
        self._write('}\n')


if __name__ == '__main__':
    import sys
    print('graph:')

    def _add_vertex(id, size, attribute_names, attribute_values):
        print('vertex:', id, size, attribute_names, attribute_values)

    def _add_edge(src, dest, attribute_names, attribute_values):
        print('edge:', src, dest, attribute_names, attribute_values)

    with open('parser-examples/graph.gxt') as f:
        parse(f, _add_vertex, _add_edge)

    def _add_minhash(id, hashes):
        print('minhashes:', id, hashes)

    print()
    print('minhashes:')

    with open('parser-examples/graph.gxt.mxt') as f:
        parse_minhash(f, _add_minhash)

    print()
    print('output gxt:')

    f = sys.stdout
    writer = Writer(f, ['label', 'fill'], ['label'])
    writer.add_vertex(1, 3, ['foo', 'red'])
    writer.add_vertex(2, 2, ['bar', 'blue'])
    writer.add_vertex(3, 1, ['baz', 'black'])
    writer.add_edge(1, 2, ['a'])
    writer.add_edge(2, 3, ['b'])
    writer.done()

    print()
    print('output gml:')

    writer = GmlWriter(f, ['label', 'fill'], ['label'])
    writer.add_vertex(1, 3, ['foo', 'red'])
    writer.add_vertex(2, 2, ['bar', 'blue'])
    writer.add_vertex(3, 1, ['baz', 'black'])
    writer.add_edge(1, 2, ['a'])
    writer.add_edge(2, 3, ['b'])
    writer.done()

    print()
    print('output dot:')

    writer = DotWriter(f)
    writer.add_vertex(1)
    writer.add_vertex(2)
    writer.add_vertex(3)
    writer.add_edge(1, 2)
    writer.add_edge(2, 3)
    writer.done()
