#! /usr/bin/env python3

import unittest

from spacegraphcats.graph import Graph
from spacegraphcats.rdomset import ldo, rdomset, test_domset
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

		# Test dominating set
		domset, _ = rdomset(g,1)

		# Test whether domset is indeed a domset
		for x in g:
			dominated = x in domset
			for y in g.neighbours(x):
				dominated = dominated or y in domset
			self.assertTrue(dominated)

		# Let's make sure the debug function in rdomset agrees
		self.assertTrue(test_domset(g, domset, 1))

		# Test higher-radius dominating sets
		for r in range(2,5):
			domset, _ = rdomset(g, r)
			self.assertTrue(test_domset(g, domset, r))

if __name__ == '__main__':
    unittest.main()
