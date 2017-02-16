#!/usr/bin/env python3

class IdentityHash:
    def __init__(self):
        pass
    def __getitem__(self, item):
        return item

"""Parser for simple graph format."""


def _parse_line(line):
    return list(map(str.strip, line.split(',')))


def _mapstr(items):
    return list(map(str, items))


def parse(graph_file, add_vertex=None, add_edge=None, consecutive_ids=False):
    """Parse a graph and call provided methods with vertices and edges."""
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
            add_vertex(_get_consecutive_id(int(parsed[0])), int(parsed[1]),
                       vertex_attributes, parsed[2:])
        next_line = graph_file.readline()

    if add_edge is None:
        # we won't be doing anything with the edges anyway
        return id_map;

    # read edges
    edge_attributes = _parse_line(graph_file.readline())[2:]

    next_line = graph_file.readline()
    while len(next_line) > 1:
        parsed = _parse_line(next_line)
        add_edge(_get_consecutive_id(int(parsed[0])), _get_consecutive_id(int(parsed[1])),
                 edge_attributes, parsed[2:])
        next_line = graph_file.readline()

    return id_map


def parse_minhash(minhash_file, add_minhash):
    """Parse minhash (.mxt) file."""
    for line in minhash_file:
        if len(line) < 2:
            continue
        parsed = _parse_line(line)
        add_minhash(parsed[0], list(map(int, map(str.strip, parsed[1].split()))))

def _parse_edgelist(graph_file, add_edge):
    """Parse and edgelist (.ext) file."""
    for line in graph_file:
        parsed = _parse_line(line)
        add_edge(parsed[0], parsed[1])

def write_edgelist(graph_file, edges):
    """Write an edgelist into an .ext file."""
    for u, v in edges:
        graph_file.write('{},{}\n'.format(u, v))

class Writer(object):
    """Writer for the gxt graph format.

    You need to either pass the vertex and edge attributes (the names)
    when you initialize the writer or when you add a vertex or edge.
    """

    def __init__(self, graph_file, vertex_attributes=None, edge_attributes=None):
        """Initialize graph writer."""
        self.file = graph_file
        self.vertex_header_written = False
        self.edge_header_written = False

        if vertex_attributes is not None:
            self.vertex_header = ','.join(['id', 'size'] + vertex_attributes)

        if edge_attributes is not None:
            self.edge_header = ','.join(['src', 'dest'] + edge_attributes)

    def add_vertex(self, vertex_id, size, attribute_values=None, vertex_attributes=None):
        """Add a vertex to the output. Don't add edges after adding nodes."""
        if attribute_values is None:
            attribute_values = []

        assert not self.edge_header_written

        if not hasattr(self, 'vertex_header'):
            self.vertex_header = ','.join(['id', 'size'] + vertex_attributes)

        if not self.vertex_header_written:
            self.file.write(self.vertex_header + '\n')
            self.vertex_header_written = True
        self.file.write(','.join(_mapstr([vertex_id, size] + attribute_values)) + '\n')

    def add_edge(self, src, dest, attribute_values=None, edge_attributes=None):
        """Add an edge to the output. Add all the nodes before you add edges."""
        if attribute_values is None:
            attribute_values = []

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


class GmlWriter(object):
    """Similar to the writer for gxt above but for gml."""

    def __init__(self, graph_file, vertex_attributes=None, edge_attributes=None, directed=False):
        """Initialize graph writer."""
        self.file = graph_file

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
        try:
            return int(value)
        except ValueError:
            pass
        return '"{}"'.format(value)

    def add_vertex(self, vertex_id, size, attribute_values=None, vertex_attributes=None):
        """Add a vertex to the output."""
        if attribute_values is None:
            attribute_values = []

        if not hasattr(self, 'vertex_attributes'):
            self.vertex_attributes = vertex_attributes

        self._write('  node [\n')
        self._write('    id {}\n'.format(vertex_id))
        self._write('    size {}\n'.format(size))
        for k, v in zip(self.vertex_attributes, attribute_values):
            self._write('    {} {}\n'.format(k, self._quote(v)))
        self._write('  ]\n')

    def add_edge(self, src, dest, attribute_values=None, edge_attributes=None):
        """Add an edge to the output."""
        if attribute_values is None:
            attribute_values = []

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


class DotWriter(object):
    """Similar to the writer for gxt above but for dot."""

    def __init__(self, graph_file):
        """Initialize graph writer."""
        self.file = graph_file

        self._write('graph G {\n')

    def _write(self, string):
        self.file.write(string)

    def add_vertex(self, vertex_id):
        """Add a vertex to the output."""
        self._write('  {};\n'.format(vertex_id))

    def add_edge(self, src, dest):
        """Add an edge to the output."""
        self._write('  {} -- {};\n'.format(src, dest))

    def done(self):
        """Call when done."""
        self._write('}\n')


if __name__ == '__main__':
    import sys
    print('graph:')

    def _add_vertex(vertex_id, size, attribute_names, attribute_values):
        print('vertex:', vertex_id, size, attribute_names, attribute_values)

    def _add_edge(src, dest, attribute_names, attribute_values):
        print('edge:', src, dest, attribute_names, attribute_values)

    with open('parser-examples/graph.gxt') as f:
        parse(f, _add_vertex, _add_edge)

    def _add_minhash(vertex_id, hashes):
        print('minhashes:', vertex_id, hashes)

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

    gml_writer = GmlWriter(f, ['label', 'fill'], ['label'])
    gml_writer.add_vertex(1, 3, ['foo', 'red'])
    gml_writer.add_vertex(2, 2, ['bar', 'blue'])
    gml_writer.add_vertex(3, 1, ['baz', 'black'])
    gml_writer.add_edge(1, 2, ['a'])
    gml_writer.add_edge(2, 3, ['b'])
    gml_writer.done()

    print()
    print('output dot:')

    dot_writer = DotWriter(f)
    dot_writer.add_vertex(1)
    dot_writer.add_vertex(2)
    dot_writer.add_vertex(3)
    dot_writer.add_edge(1, 2)
    dot_writer.add_edge(2, 3)
    dot_writer.done()
