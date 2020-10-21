#! /usr/bin/env python
"""
Benchmark the rdomset (catlas level 1) algorithm, without I/O considerations.
"""
import sys
import os

# add spacegraphcats package to import path:
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import spacegraphcats
from spacegraphcats.catlas import catlas
import argparse
import sys
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("project", help="Project directory", type=str)
    parser.add_argument("radius", help="Catlas radius", type=int)
    parser.add_argument("-o", "--output", type=str)
    args = parser.parse_args()

    proj = catlas.Project(args.project, args.radius, False)

    # load the graph
    proj.load_furthest_checkpoint()
    nodes_in_layer_0 = len(proj.graph)

    # build the first layer only, for rdomset benchmarking.
    start = time.time()
    catlas.CAtlas.build(proj, benchmark_only=True)
    end = time.time()

    outfp = sys.stdout
    if args.output:
        outfp = open(args.output, "at")
    print(
        "{},{},{:.1f},{},rdomset".format(
            nodes_in_layer_0, args.radius, end - start, args.project
        ),
        file=outfp,
    )


if __name__ == "__main__":
    sys.exit(main())
