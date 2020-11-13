import unittest
import itertools
import random

from spacegraphcats.catlas.graph import Graph
from spacegraphcats.catlas.rdomset import low_degree_orientation


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


if __name__ == "__main__":
    unittest.main()
