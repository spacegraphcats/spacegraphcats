import unittest
import itertools, random        

from spacegraphcats.graph import Graph
from spacegraphcats.rdomset import low_degree_orientation

class ParserRDomset(unittest.TestCase):
    def test_ldo(self):
        # Test on complete graph
        n = 25
        g = Graph(num_nodes=n, radius=3)
        for x,y in itertools.combinations(range(n), 2):
            g.add_arc(x,y)
            g.add_arc(y,x)
        
        low_degree_orientation(g)

        for x,y in itertools.combinations(range(n), 2):
            # Every edge must be present as an arc
            self.assertTrue( g.adjacent(x,y) or g.adjacent(y,x) )

        # The low-degree orientation of a simple graph must be simple 
        # and every arc has weight 1. 
        self.assertEqual(g.arcs(2), [])
        self.assertEqual(len(g.arcs(1)), 300)

        # Test on (sparse) random graph
        n = 100
        d = 5
        p = d / n

        g = Graph(num_nodes=n, radius=3)
        count = 0
        for x, y in itertools.combinations(range(n), 2):
            if random.random() < p:
                if (x != y):
                    count += 1
                    g.add_arc(x,y)
                    g.add_arc(y,x)

        low_degree_orientation(g)

        self.assertEqual(g.arcs(2), [])
        self.assertEqual(len(g.arcs(1)), count)


if __name__ == '__main__':
    unittest.main()
