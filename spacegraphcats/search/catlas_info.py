#! /usr/bin/env python
import argparse
import os
import sys

from .search_utils import (load_dag, load_layer1_to_cdbg)
from ..catlas.graph_io import read_from_gxt

def main(argv=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix')
    p.add_argument('-k', '--ksize', default=21, type=int,
                   help='list of k-mer sizes (default: 31)')
    args = p.parse_args(argv)

    assert args.catlas_prefix.split('_')[-1] == 'r1'

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')
    contigfile = os.path.join(args.catlas_prefix, "contigs.fa.gz")
    gxtfile = os.path.join(args.catlas_prefix, 'cdbg.gxt')

    top_node_id, dag, dag_up, dag_levels, cdbg_to_catlas = load_dag(catlas)
    print('loaded {} nodes and {} layers from catlas {}'.format(len(dag), dag_levels[top_node_id], catlas))

    print('top catlas node {} has {} children.'.format(top_node_id,
                                                       len(dag[top_node_id])))

    layer1_to_cdbg = load_layer1_to_cdbg(cdbg_to_catlas, domfile)
    x = set()
    for v in layer1_to_cdbg.values():
        x.update(v)
    total_cdbg_count = len(x)
    print('{} layer 1 catlas nodes, corresponding to {} cDBG nodes.'.format(len(layer1_to_cdbg), len(x)))

    with open(gxtfile, 'rt') as fp:
        graph = read_from_gxt(fp, 1, False)

    num_arcs = graph.num_arcs()

    print('{} nodes, {} arcs, {:.1f} average.'.format(len(graph), num_arcs, num_arcs / len(graph)))

    # @CTB: could put size distribution of those nodes here...?


if __name__ == '__main__':
    # import cProfile
    # cProfile.run('main()', 'search_stats')

    main()
