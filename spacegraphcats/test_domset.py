#! /usr/bin/env python3

import unittest

from spacegraphcats.graph import Graph

class DomsetTest(unittest.TestCase):

	def test_domset(self):
		# g = Graph()
		# g.add_node(0)
		# g.add_edge(1,2)
		# g.add_edge(3,4)
		# g.add_edge(5,6).add_edge(6,7).add_edge(7,5)
		# g.add_edge(8,9).add_edge(8,10).add_edge(8,11)

		# self.assertTrue(g.num_components() == 5)

		# comps = g.components()
		# self.assertTrue(len(comps) == 5)

		# union = set()
		# for c in comps:
		# 	union |= set(c)
		self.assertTrue(len(union) == len(g))

	# def test_ldo(self):
	# 	from spacegraphcats.rdomset import ldo
	# 	import itertools, random
		
	# 	# Test on complete graph
	# 	n = 25
	# 	g = Graph.on(list(range(n)))
	# 	for x,y in itertools.combinations(range(n),2):
	# 		g.add_edge(x,y)
	# 	tfgraph = ldo(g)

	# 	for x,y in itertools.combinations(range(n),2):
	# 		# Every edge must be present as an arc
	# 		self.assertTrue( tfgraph.adjacent(x,y) or tfgraph.adjacent(y,x) )
	# 		# The low-degree orientation of a simple graph must be simple 
	# 		# and every arc has weight 1. 
	# 		self.assertTrue( tfgraph.weight(x,y) == 1 or tfgraph.weight(y,x) == 1 )
	# 		self.assertTrue( tfgraph.weight(x,y) == None or tfgraph.weight(y,x) == None )

	# 	# Test on (sparse) random graph
	# 	n = 100
	# 	d = 5
	# 	p = d / n
	# 	print(p)

	# 	g = Graph.on(list(range(n)))
	# 	for x,y in itertools.combinations(range(n),2):
	# 		if random.random() < p:
	# 			g.add_edge(x,y)
	# 	g.remove_loops()
	# 	tfgraph = ldo(g)

	# 	for x,y in g.edges():
	# 		# Every edge must be present as an arc
	# 		self.assertTrue( tfgraph.adjacent(x,y) or tfgraph.adjacent(y,x) )
	# 		# The low-degree orientation of a simple graph must be simple 
	# 		# and every arc has weight 1. 
	# 		self.assertTrue( tfgraph.weight(x,y) == 1 or tfgraph.weight(y,x) == 1 )
	# 		self.assertTrue( tfgraph.weight(x,y) == None or tfgraph.weight(y,x) == None )			

if __name__ == '__main__':
    unittest.main()
