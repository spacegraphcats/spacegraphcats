#! /usr/bin/env python3

import unittest
import itertools
from spacegraphcats.catlas.graph import Graph


class GraphTest(unittest.TestCase):
    def test_graph(self):
        graph = Graph(num_nodes=5, radius=5)

        # Query in-neighbourhood by weight
        graph.add_arc(1, 0, 1).add_arc(2, 0, 1).add_arc(3, 0, 1).add_arc(4, 0, 1)
        neighbors = set(graph.in_neighbors(0, 1))
        self.assertEqual(neighbors, set([1, 2, 3, 4]))

        neighbors = set(graph.in_neighbors(0))
        self.assertEqual(neighbors, set([(1, 1), (2, 1), (3, 1), (4, 1)]))

        self.assertEqual(graph.num_arcs(), 4)
        self.assertEqual(graph.num_arcs(1), 4)
        self.assertEqual(graph.num_arcs(2), 0)

        self.assertEqual(graph.in_degree(0), 4)
        self.assertEqual(graph.in_degree(1), 0)

    def test_frat_pairs(self):
        frat_graph = Graph(num_nodes=9, radius=5)

        frat_graph.add_arc(1, 0, 1)
        frat_graph.add_arc(2, 0, 1)
        frat_graph.add_arc(3, 0, 1)
        frat_graph.add_arc(4, 0, 1)
        frat_graph.add_arc(5, 0, 5)
        frat_graph.add_arc(6, 0, 5)
        frat_graph.add_arc(7, 0, 5)
        frat_graph.add_arc(8, 0, 5)
        # need to make sure we don't include adjacent pairs
        frat_graph.add_arc(8, 1, 1)

        # arcs of weight 1
        w1 = [i for i in range(1, 5)]
        # arcs of weight 5
        w5 = [i for i in range(5, 9)]
        # frat pairs of weight 2
        frat_pairs2 = set(itertools.combinations(w1, 2))
        # frat pairs of weight 6
        frat_pairs6 = set(itertools.product(w1, w5))
        frat_pairs6.remove((1, 8))

        self.assertEqual(set(frat_graph.fraternal_pairs(0, 2)), frat_pairs2)
        self.assertEqual(set(frat_graph.fraternal_pairs(0, 6)), frat_pairs6)

    def test_trans_pairs(self):
        graph = Graph(num_nodes=9, radius=3)

        graph.add_arc(1, 0, 2)
        graph.add_arc(2, 0, 3)
        graph.add_arc(3, 0, 1)
        graph.add_arc(4, 1, 2)
        graph.add_arc(5, 1, 2)
        graph.add_arc(6, 0, 1)
        graph.add_arc(6, 2, 1)
        graph.add_arc(7, 3, 3)
        graph.add_arc(8, 3, 2)

        trans_pairs = set([(7, 0), (4, 0), (5, 0)])
        self.assertEqual(set(graph.transitive_pairs(0, 4)), trans_pairs)


if __name__ == "__main__":
    unittest.main()
