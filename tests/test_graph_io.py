import unittest
from io import StringIO

from spacegraphcats.catlas.graph_io import read_from_gxt, write_to_gxt
from spacegraphcats.catlas.graph import Graph


class IOTest(unittest.TestCase):
    def test_writing_and_reading(self):
        f = StringIO()

        graph = Graph(5)
        graph.add_arc(1, 0).add_arc(2, 0).add_arc(3, 0).add_arc(4, 0)

        write_to_gxt(f, graph, 1)
        f.seek(0)

        parsed = read_from_gxt(f, 5, True)

        self.assertEqual(list(parsed.arcs()), list(graph.arcs()))
        self.assertEqual(len(parsed), len(graph))

    def test_writing_and_reading_no_weight(self):
        f = StringIO()

        graph = Graph(5)
        graph.add_arc(1, 0).add_arc(2, 0).add_arc(3, 0).add_arc(4, 0)

        write_to_gxt(f, graph)
        f.seek(0)

        parsed = read_from_gxt(f, 5, True)

        self.assertEqual(list(parsed.arcs()), list(graph.arcs()))
        self.assertEqual(len(parsed), len(graph))


if __name__ == "__main__":
    unittest.main()
