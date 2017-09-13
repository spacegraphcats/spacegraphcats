#! /usr/bin/env python
import sys
import os
import argparse
import screed
import collections
from pickle import load, dump
from search import search_utils
from search_utils import (load_dag, load_layer0_to_cdbg)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('cover')
    p.add_argument('reads_inp')
    p.add_argument('reads_outp')
    args = p.parse_args()

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # load catlas DAG
    top_node_id, dag, dag_up, dag_levels = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    # load mapping between dom nodes and cDBG/graph nodes:
    layer0_to_cdbg = load_layer0_to_cdbg(catlas, domfile)
    print('loaded {} layer 0 catlas nodes'.format(len(layer0_to_cdbg)))
    x = set()
    for v in layer0_to_cdbg.values():
        x.update(v)
    print('...corresponding to {} cDBG nodes.'.format(len(x)))

    # reverse
    cdbg_to_layer0 = collections.defaultdict(set)
    for k, v in layer0_to_cdbg.items():
        for vv in v:
            cdbg_to_layer0[vv].add(k)

    covernodes = load(open(args.cover, 'rb'))

    def get_leaf_nodes(catlas_id):
        x = set()
        if dag[catlas_id]:
            for child in dag[catlas_id]:
                x.update(get_leaf_nodes(child))
        else:
            x.add(catlas_id)

        return x

    cdbg_to_cover = collections.defaultdict(set)
    total_cdbg_covered = set()
    for cover_node in covernodes:
        # get all leaf nodes under given cover node
        leaf_nodes = get_leaf_nodes(cover_node)
        cdbg_nodes = set()

        # convert all leaf nodes to cDBG IDs.
        for l in leaf_nodes:
            cdbg_nodes.update(cdbg_to_layer0.get(l))

        # now construct reverse index: cDBG node to cover_node
        for c in cdbg_nodes:
            cdbg_to_cover[c].add(cover_node)

        total_cdbg_covered.update(cdbg_nodes)

    print(len(cdbg_to_cover), len(total_cdbg_covered))

    total_bp = 0
    watermark_size = 1e7
    watermark = watermark_size
    no_tags = 0
    no_tags_bp = 0
    total_seqs = 0

    n = 0
    outfp = open(args.reads_outp, 'wt')
    for record in screed.open(args.reads_inp):
        n += 1
        if total_bp >= watermark:
            print('... {:5.2e} bp thru reads'.format(int(watermark)),
                  file=sys.stderr)
            watermark += watermark_size

        total_bp += len(record.sequence)
        total_seqs += 1

        if not record.name:
            continue

        labels = map(int, record.name.split(','))

        for t in labels:
            if cdbg_to_cover.get(t):
                outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
                break

    print(no_tags, total_seqs)
    print('{:5.2e} {:5.2e}'.format(no_tags_bp, total_bp))


if __name__ == '__main__':
    main()
