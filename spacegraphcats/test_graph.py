#! /usr/bin/env python3

import unittest

from spacegraphcats.graph import Graph, TFGraph

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

	def test_tfgraph(self):
		tfgraph = TFGraph(list(range(5)))

		# Query in-neighbourhood by weight
		tfgraph.add_arc(1,0,1).add_arc(2,0,1).add_arc(3,0,1).add_arc(4,0,1)
		N = set(tfgraph.in_neighbours_weight(0,1))
		self.assertTrue(N == set([1,2,3,4]))

		# There should be no transitive triples in the current graph:
		# All arcs point to vertex 0 which has not outgoing arc
		for x in tfgraph:
			trips = list(tfgraph.trans_trips(x))
			self.assertTrue(len(trips) == 0)

		# No vertex except 0 should produce fraternal triples
		for x in [1,2,3,4]:
			trips = list(tfgraph.frat_trips(x))
			self.assertTrue(len(trips) == 0)			

		# Every pair of [1,2,3,4] is fraternal via vertex 0 with
		# weight 2
		trips = list(tfgraph.frat_trips(0))
		self.assertTrue(len(trips) == 6)
		for _,_,w in trips:
			self.assertTrue(w == 2)

		# Add an arc of weight 5 from 0 to 1. Now the vertices 
		# [1,2,3,4] are in the second back-neighbourhood of 1 with 
		# total weight 6
		tfgraph.add_arc(0,1,5)
		trips = list(tfgraph.trans_trips(1))
		self.assertTrue(len(trips) == 4)
		deplete = set([1,2,3,4]) 
		for x,y,w in trips:
			self.assertTrue(x in deplete)
			deplete.remove(x)
			self.assertTrue(y == 1)
			self.assertTrue(w == 6)
		self.assertTrue(len(deplete) == 0)


	def test_ldo(self):
		from spacegraphcats.rdomset import ldo
		import itertools, random
		
		# Test on complete graph
		n = 25
		g = Graph.on(list(range(n)))
		for x,y in itertools.combinations(range(n),2):
			g.add_edge(x,y)
		tfgraph = ldo(g)

		for x,y in itertools.combinations(range(n),2):
			# Every edge must be present as an arc
			self.assertTrue( tfgraph.adjacent(x,y) or tfgraph.adjacent(y,x) )
			# The low-degree orientation of a simple graph must be simple 
			# and every arc has weight 1. 
			self.assertTrue( tfgraph.weight(x,y) == 1 or tfgraph.weight(y,x) == 1 )
			self.assertTrue( tfgraph.weight(x,y) == None or tfgraph.weight(y,x) == None )

		# Test on (sparse) random graph
		n = 100
		d = 5
		p = d / n

		g = Graph.on(list(range(n)))
		for x,y in itertools.combinations(range(n),2):
			if random.random() < p:
				g.add_edge(x,y)
		g.remove_loops()
		tfgraph = ldo(g)

		for x,y in g.edges():
			# Every edge must be present as an arc
			self.assertTrue( tfgraph.adjacent(x,y) or tfgraph.adjacent(y,x) )
			# The low-degree orientation of a simple graph must be simple 
			# and every arc has weight 1. 
			self.assertTrue( tfgraph.weight(x,y) == 1 or tfgraph.weight(y,x) == 1 )
			self.assertTrue( tfgraph.weight(x,y) == None or tfgraph.weight(y,x) == None )			

if __name__ == '__main__':
    unittest.main()
