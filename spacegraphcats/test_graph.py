#! /usr/bin/env python3

import unittest

from spacegraphcats.graph import Graph

class GraphTest(unittest.TestCase):

	def test_components(self):
		g = Graph()
		g.add_node(0)
		g.add_edge(1,2)
		g.add_edge(3,4)
		g.add_edge(5,6).add_edge(6,7).add_edge(7,5)
		g.add_edge(8,9).add_edge(8,10).add_edge(8,11)

		self.assertTrue(g.num_components() == 5)

		comps = g.components()
		self.assertTrue(len(comps) == 5)

		union = set()
		for c in comps:
			union |= set(c)

		self.assertTrue(len(union) == len(g))

if __name__ == '__main__':
    unittest.main()
