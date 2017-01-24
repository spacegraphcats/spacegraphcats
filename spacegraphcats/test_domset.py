#! /usr/bin/env python3

import unittest

from spacegraphcats.graph import Graph
from spacegraphcats.rdomset import ldo, better_dvorak_reidl
import itertools, random

class DomsetTest(unittest.TestCase):

	def test_domset(self):
		n = 100
		d = 5
		p = d / n
		g = Graph.on(list(range(n)))

		for x,y in itertools.combinations(range(n),2):
			if random.random() < p:
				g.add_edge(x,y)
		g.remove_loops()
		tfgraph = ldo(g)

		domset = better_dvorak_reidl(tfgraph,1)
		# Test whether domset is indeed a domset
		for x in g:
			dominated = x in domset
			for y in g.neighbours(x):
				dominated = dominated or y in domset
			self.assertTrue(dominated)




if __name__ == '__main__':
    unittest.main()
