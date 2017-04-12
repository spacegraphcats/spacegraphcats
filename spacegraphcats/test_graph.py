#! /usr/bin/env python3

import unittest

from spacegraphcats.graph import Graph


class GraphTest(unittest.TestCase):

    def test_graph(self):
        graph = Graph(num_nodes=5, radius=5)

        # Query in-neighbourhood by weight
        graph.add_arc(1, 0, 1).add_arc(2, 0, 1).add_arc(3, 0, 1).add_arc(4, 0,
                                                                         1)
        neighbors = set(graph.in_neighbors(0, 1))
        self.assertEquals(neighbors, set([1, 2, 3, 4]))

        neighbors = set(graph.in_neighbors(0))
        self.assertEquals(neighbors, set([(1, 1), (2, 1), (3, 1), (4, 1)]))

        self.assertEquals(graph.num_arcs(), 4)
        self.assertEquals(graph.num_arcs(1), 4)
        self.assertEquals(graph.num_arcs(2), 0)

        self.assertEquals(graph.in_degree(0), 4)
        self.assertEquals(graph.in_degree(1), 0)


if __name__ == '__main__':
    unittest.main()
