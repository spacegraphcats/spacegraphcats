#! /usr/bin/env python3

import unittest
import os

import spacegraphcats.graph_parser as parser

DIR = os.path.dirname(os.path.realpath(__file__))


class ParserTest(unittest.TestCase):
    def test_parsing(self):
        all_ids = []
        all_sizes = []
        all_names_v = []
        all_values_v = []

        def collect_vertex(id, size, names, values):
            all_ids.append(id)
            all_sizes.append(size)
            all_names_v.append(names)
            all_values_v.append(values)

        all_srcs = []
        all_dests = []
        all_names_e = []
        all_values_e = []

        def collect_edge(src, dest, names, values):
            all_srcs.append(src)
            all_dests.append(dest)
            all_names_e.append(names)
            all_values_e.append(values)

        with open(os.path.join(DIR, 'parser-examples/graph.gxt')) as f:
            self.p = parser.parse(f, collect_vertex, collect_edge)

        self.assertEqual(all_ids, [1, 2, 3])
        self.assertEqual(all_sizes, [2, 3, 1])
        self.assertEqual(all_names_v, [['label', 'fill'], ['label', 'fill'], ['label', 'fill']])
        self.assertEqual(all_values_v, [['foo', 'red'], ['bar', 'green'], ['batman', 'black']])

        self.assertEqual(all_srcs, [1, 2])
        self.assertEqual(all_dests, [2, 3])
        self.assertEqual(all_names_e, [['label'], ['label']])
        self.assertEqual(all_values_e, [['a'], ['b']])

if __name__ == '__main__':
    unittest.main()
