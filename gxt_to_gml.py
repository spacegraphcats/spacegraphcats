#!/usr/bin/env python3

import sys
import argparse

from spacegraphcats.graph_parser import parse, GmlWriter

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='gxt input file', nargs='?',
                        type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('output', help='gml output file', nargs='?',
                        type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--directed', help='set flag if graph is supposed to be directed', action='store_true', default=False)
    args = parser.parse_args()

    writer = GmlWriter(args.output, directed=args.directed)

    def add_vertex(v, size, attributes, values):
        writer.add_vertex(v, size, values, attributes)

    def add_edge(u, v, attributes, values):
        writer.add_edge(u, v, values, attributes)

    parse(args.input, add_vertex, add_edge)
    writer.done()
