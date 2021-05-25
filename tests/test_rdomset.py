import unittest
import itertools
import random

from spacegraphcats.catlas.graph import Graph
from spacegraphcats.catlas.rdomset import low_degree_orientation, domination_graph
from sortedcontainers import SortedSet, SortedDict


class ParserRDomset(unittest.TestCase):
    def test_ldo_random(self):
        # Test on complete graph
        n = 25
        g = Graph(num_nodes=n, radius=3)
        for x, y in itertools.combinations(range(n), 2):
            g.add_arc(x, y)
            g.add_arc(y, x)

        low_degree_orientation(g)

        for x, y in itertools.combinations(range(n), 2):
            # Every edge must be present as an arc
            self.assertTrue(g.adjacent(x, y) or g.adjacent(y, x))

        # The low-degree orientation of a simple graph must be simple
        # and every arc has weight 1.
        self.assertEqual(len(list(g.arcs(2))), 0)
        self.assertEqual(len(list(g.arcs(1))), 300)

        # Test on (sparse) random graph
        n = 100
        d = 5
        p = d / n

        g = Graph(num_nodes=n, radius=3)
        count = 0
        for x, y in itertools.combinations(range(n), 2):
            if random.random() < p:
                if x != y:
                    count += 1
                    g.add_arc(x, y)
                    g.add_arc(y, x)

        low_degree_orientation(g)

        self.assertEqual(len(list(g.arcs(2))), 0)
        self.assertEqual(len(list(g.arcs(1))), count)

    def test_ldo_example(self):
        g = Graph(num_nodes=10, radius=3)

        edges = [
            (0, 6),
            (6, 7),
            (6, 5),
            (6, 8),
            (5, 8),
            (8, 9),
            (8, 2),
            (2, 1),
            (2, 3),
            (3, 4),
        ]

        for x, y in edges:
            g.add_arc(x, y)
            g.add_arc(y, x)

        low_degree_orientation(g)

        self.assertEqual(len(list(g.arcs(2))), 0)
        self.assertEqual(len(list(g.arcs(1))), 10)

        # unambiguous arcs
        self.assertTrue(
            set([(6, 0), (2, 1), (8, 2), (2, 3), (3, 4), (6, 7), (8, 9)]).issubset(
                set(g.arcs(1))
            )
        )

    def test_domination_graph(self):
        r = 3
        g = Graph(num_nodes=10, radius=r)

        edges = [(0, 6), (6, 7), (6, 5), (6, 8), (5, 8), (8, 9), (8, 2), (2, 1), (2, 3), (3, 4)]

        for x, y in edges:
            g.add_arc(x, y)
            g.add_arc(y, x)

        low_degree_orientation(g)

        # normal cases

        # Vertices must be assigned to their closest dominators.
        # There are no ties with the following r-dominating sets.
        domgraph, dominated = domination_graph(g, [0, 4], r)
        self.assertEqual(
            dominated,
            SortedDict({0: SortedSet([0, 5, 6, 7, 8, 9]), 4: SortedSet([1, 2, 3, 4])})
        )
        self.assertEqual(set(domgraph.arcs(1)), set([(0, 4), (4, 0)]))

        domgraph, dominated = domination_graph(g, [3, 6], r)
        self.assertEqual(
            dominated,
            SortedDict({3: SortedSet([1, 2, 3, 4]), 6: SortedSet([0, 5, 6, 7, 8, 9])})
        )
        self.assertEqual(set(domgraph.arcs(1)), set([(3, 6), (6, 3)]))

        # Vertices 8, 9 are equidistance from the dominators, but they should be assigned to
        # the earliest dominator, which is vertex 2. Thus, the resulting pieces should induce
        # connected subgraphs.
        domgraph, dominated = domination_graph(g, [2, 5, 6], r)
        self.assertEqual(
            dominated,
            SortedDict({2: SortedSet([1, 2, 3, 4, 8, 9]), 5: SortedSet([5]), 6: SortedSet([0, 6, 7])})
        )
        self.assertEqual(set(domgraph.arcs(1)), set([(2, 5), (2, 6), (5, 2), (5, 6), (6, 2), (6, 5)]))

        # error cases (should raise an error when the given vertex set is not an r-dominating set)
        self.assertRaises(AssertionError, lambda: domination_graph(g, [], r))
        self.assertRaises(AssertionError, lambda: domination_graph(g, [0], r))
        self.assertRaises(AssertionError, lambda: domination_graph(g, [0, 5, 6, 7, 9], r))


if __name__ == "__main__":
    unittest.main()
