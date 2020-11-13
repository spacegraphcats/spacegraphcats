#! /usr/bin/env python
import argparse
import os
import sys

from .catlas import CAtlas
from ..catlas.graph_io import read_from_gxt


def main(argv=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument("catlas_prefix")
    args = p.parse_args(argv)

    assert args.catlas_prefix.split("_")[-1] == "r1"

    basename = os.path.basename(args.catlas_prefix)
    gxtfile = os.path.join(args.catlas_prefix, "cdbg.gxt")

    catlas = CAtlas(args.catlas_prefix)
    top_node_id, dag, dag_levels = catlas.root, catlas.children, catlas.levels

    print(
        "loaded {} nodes and {} layers from catlas {}".format(
            len(dag), dag_levels[top_node_id], args.catlas_prefix
        )
    )

    print(
        "top catlas node {} has {} children.".format(top_node_id, len(dag[top_node_id]))
    )

    layer1_to_cdbg = catlas.layer1_to_cdbg
    x = set()
    for v in layer1_to_cdbg.values():
        x.update(v)
    total_cdbg_count = len(x)
    print(
        "{} layer 1 catlas nodes, corresponding to {} cDBG nodes.".format(
            len(layer1_to_cdbg), len(x)
        )
    )

    with open(gxtfile, "rt") as fp:
        graph = read_from_gxt(fp, 1, False)

    num_arcs = graph.num_arcs()

    print(
        "{} nodes, {} arcs, {:.1f} average.".format(
            len(graph), num_arcs, num_arcs / len(graph)
        )
    )

    # CTB: could put size distribution of those nodes here...?

    return 0


if __name__ == "__main__":
    # import cProfile
    # cProfile.run('main()', 'search_stats')

    sys.exit(main())
