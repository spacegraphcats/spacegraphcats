import unittest
import os
from io import StringIO

from spacegraphcats.catlas import graph_parser as parser

DIR = os.path.dirname(os.path.realpath(__file__))


class ParserTest(unittest.TestCase):
    def test_graph_parsing(self):
        num_vertices = []

        def create_vertices(n):
            num_vertices.append(n)

        all_edges = []

        def collect_edge(src, dest):
            all_edges.append((src, dest))

        with open(os.path.join(DIR, "test-data/parser-examples/graph.gxt")) as f:
            id_map = parser.parse(f, create_vertices, collect_edge)

        self.assertEqual(num_vertices, [3])
        self.assertEqual(all_edges, ([(0, 1), (1, 2)]))

    def test_graph_writing(self):
        output = StringIO()

        parser.write(output, 3, [(0, 1), (1, 2)])

        self.assertEqual(output.getvalue(), "3\n0 1\n1 2\n")

    def test_mxt_parsing(self):
        minhashes = {}

        def collect_minhashes(vertex_id, hashes):
            minhashes[vertex_id] = hashes

        with open(os.path.join(DIR, "test-data/parser-examples/graph.gxt.mxt")) as f:
            parser.parse_minhash(f, collect_minhashes)

        self.assertEqual(
            minhashes,
            {
                0: [14891351629450783567, 8602412685556304666, 15005322196398795210],
                2: [17662871537941316484],
            },
        )


if __name__ == "__main__":
    unittest.main()
