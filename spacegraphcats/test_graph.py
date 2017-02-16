#! /usr/bin/env python3

import unittest

from spacegraphcats.graph import Graph

class GraphTest(unittest.TestCase):

    def test_graph(self):
        graph = Graph(num_nodes=5, radius=5)

        # Query in-neighbourhood by weight
        graph.add_arc(1, 0, 1).add_arc(2, 0, 1).add_arc(3, 0, 1).add_arc(4, 0, 1)
        neighbors = set(graph.in_neighbors(0, 1))
        self.assertEquals(neighbors, set([1, 2, 3, 4]))

        neighbors = set(graph.in_neighbors(0))
        self.assertEquals(neighbors, set([(1, 1), (2, 1), (3, 1), (4, 1)]))

        self.assertEquals(graph.num_arcs(), 4)
        self.assertEquals(graph.num_arcs(1), 4)
        self.assertEquals(graph.num_arcs(2), 0)

        self.assertEquals(graph.in_degree(0), 4)
        self.assertEquals(graph.in_degree(1), 0)

    def test_components(self):
        g = Graph(num_nodes=12)
        g.add_arc(1,2)
        g.add_arc(3,4)
        g.add_arc(5,6).add_arc(6,7).add_arc(7,5)
        g.add_arc(8,9).add_arc(8,10).add_arc(8,11)

        self.assertEquals(g.num_components(), 5)

        comps = g.components()
        self.assertEquals(len(comps), 5)

        union = set()
        for c in comps:
            union |= set(c)

        self.assertEquals(len(union), len(g))
        

if __name__ == '__main__':
    unittest.main()
