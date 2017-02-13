#! /usr/bin/env python3

import unittest

from spacegraphcats.graph import Graph
from spacegraphcats.rdomset import *
import itertools, random

class DomsetTest(unittest.TestCase):

    def test_rdomset(self):
        """ Testing 2-domination on a path of length four where it is pretty
            clear how the resulting data structures look like. """
        g = Graph(n=5)
        g.add_edge(0,1).add_edge(1,2).add_edge(2,3).add_edge(3,4)

        augg = dtf(g,2)

        # Verify that dominating sets are recognized properly
        self.assertTrue(verify_domset(g, [1,3], 1))
        self.assertFalse(verify_domset(g, [2], 1))
        self.assertTrue(verify_domset(g, [2], 2))
        self.assertFalse(verify_domset(g, [1], 2))

        # Verify that dominators are computed correctly
        dominators = calc_dominators(augg, [0,4], 2)
        self.assertTrue( dominators[0][0] == set([0]) ) # Vertex 0 is dominated at distance zero by itself
        self.assertTrue( dominators[0][1] == set([]) ) # Vertex 0 does not have a dominator at distance one
        self.assertTrue( dominators[0][2] == set([]) ) # Vertex 0 does not have a dominator at distance two

        self.assertTrue( dominators[1][0] == set([]) ) # Vertex 1 is not part of the dominating set
        self.assertTrue( dominators[1][1] == set([0]) ) # Vertex 1 is dominated by 0 at distance one
        self.assertTrue( dominators[1][2] == set([]) ) # Vertex 1 is not dominated at distance two

        self.assertTrue( dominators[2][0] == set([]) ) # Vertex 2 is not in the domset
        self.assertTrue( dominators[2][1] == set([]) ) # Vertex 2 has no distance-one dominators
        self.assertTrue( dominators[2][2] == set([0,4]) ) # Vertex 2 is dominated by 0 and 0 at distance two

        # The domination graph is formed on the vertices of the dominating set.
        # Two dominators are connected by an edge if they both dominate a common vertex.
        # In general this function can return a slightly larger dominating set, if the
        # resulting graph would be disconnected.
        h, assignment = calc_domination_graph(g, augg, [0,4], dominators, 2)
        self.assertTrue(h.nodes == set([0,4]))
        self.assertTrue(h.adjacent(0,4)) # Vertex 0 and 4 both dominate vertex 2
        #self.assertTrue(newdomset == set([0,4])) # There is no need to augment the domset

        # The assignment map tells us for each vertex of g what its closest dominators are.
        # We expect [0,1] to belong to 0, [3,4] to belong to 4 and only 2 to belong to both 0 and 4.
        self.assertTrue(assignment[0] == set([0]))
        self.assertTrue(assignment[1] == set([0]))
        self.assertTrue(assignment[2] == set([0,4]))
        self.assertTrue(assignment[3] == set([4]))
        self.assertTrue(assignment[4] == set([4]))

    def test_rdomset_augmentation(self):
        """ Covers a simple test case in which a dominating set must be modified
            to ensure that the resultign domination graph is connected. """
        g = Graph(n=5)
        g.add_edge(0,1).add_edge(1,2).add_edge(2,3).add_edge(3,4).add_edge(4,5)

        augg = dtf(g, 2)
        self.assertTrue(verify_domset(g, [0,5], 2)) # Verify that the endpoints of the path are a 2-domset

        # Verify that 0 and 5 are 2-scattered
        self.assertTrue(verify_scattered(g, [0,5], 2))

        # Verify this fact again by hand using dominators
        dominators = calc_dominators(augg, [0,5], 2)
        self.assertTrue(dominators[2][2] == set([0])) # Vertex 2 is only dominated by 0 at distance 2
        self.assertTrue(dominators[3][2] == set([5])) # Vertex 3 is only dominated by 5 at distance 2

        # The domgraph computation now needs to add an additional vertex in order to make
        # it connected. The new domset might not be a superset of the old domset!
        #h, assignment = calc_domination_graph(g, augg, [0,5],dominators, 2)
        #self.assertTrue( len(newdomset) > 2 )

    def test_rdomset_random(self):
        """ This test covers pretty much the whole r-domset infrastructure
            and asserts that what we compute are indeed r-dominating sets. """
        n = 100
        d = 5
        p = d / n
        g = Graph(n=n)

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
        self.assertTrue(verify_domset(g, domset, 1))

        # Test higher-radius dominating sets
        for r in range(2,5):
            domset, augg = rdomset(g, r)
            self.assertTrue(verify_domset(g, domset, r))

            # If domset is an r-domset, then the set of covered
            # vertices must be the whole graph
            covered = set(calc_dominated(augg, domset, r))
            self.assertTrue( covered == set(g.nodes) )

            # Test assignment of dominators
            dominators = calc_dominators(augg, domset, r)
            self.assertTrue(verify_dominators(g, domset, dominators, r))
            self.assertTrue(verify_dominators_strict(g, domset, dominators, r))

if __name__ == '__main__':
    unittest.main()
