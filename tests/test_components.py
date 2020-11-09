#! /usr/bin/env python3
import unittest

from spacegraphcats.catlas.graph import Graph
from spacegraphcats.catlas.components import components, num_components


class ComponentsTest(unittest.TestCase):
    def test_components(self):
        g = Graph(num_nodes=12)
        g.add_arc(1, 2)
        g.add_arc(3, 4)
        g.add_arc(5, 6).add_arc(6, 7).add_arc(7, 5)
        g.add_arc(8, 9).add_arc(8, 10).add_arc(8, 11)

        self.assertEqual(num_components(g), 5)

        comps = components(g)
        self.assertEqual(len(comps), 5)

        union = set()
        for c in comps:
            union |= set(c)

        self.assertEqual(len(union), len(g))


if __name__ == "__main__":
    unittest.main()
