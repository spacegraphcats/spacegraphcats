#! /usr/bin/env python
"""
Look for catlas nodes that have few k-mers but many cDBG nodes underneath
them, as a likely sign of strain variation.
"""
import argparse
import os
import sys
import collections
import math
import sourmash_lib
import time
import numpy

import screed
from . import search_utils
from .frontier_search import find_shadow


def main(args=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('output')
    p.add_argument('--maxsize', type=float, default=100000)
    p.add_argument('--minsize', type=float, default=5000)
    p.add_argument('-k', '--ksize', default=31, type=int,
                        help='k-mer size (default: 31)')
    args = p.parse_args(args)

    print('maxsize: {:g}'.format(args.maxsize))

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # load catlas DAG
    top_node_id, dag, dag_up, dag_levels, cdbg_to_catlas = search_utils.load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    # load mapping between dom nodes and cDBG/graph nodes:
    layer1_to_cdbg = search_utils.load_layer1_to_cdbg(cdbg_to_catlas, domfile)
    print('loaded {} layer 1 catlas nodes'.format(len(layer1_to_cdbg)))

    # find the contigs filename
    contigs = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

    # ...and kmer index.
    ki_start = time.time()
    kmer_idx = search_utils.load_kmer_index(args.catlas_prefix)
    print('loaded {} k-mers in index ({:.1f}s)'.format(len(kmer_idx.mphf_to_kmer), time.time() - ki_start))

    # calculate the k-mer sizes for each catlas node.
    node_sizes = kmer_idx.build_catlas_node_sizes(dag, dag_levels, layer1_to_cdbg)

    # find highest nodes with shadow size less than given max_size
    def find_terminal_nodes(node_id, max_size):
        node_list = set()
        for sub_id in dag[node_id]:
            # shadow size
            size = node_sizes[sub_id]

            if size < max_size:
                node_list.add(sub_id)
            else:
                children = find_terminal_nodes(sub_id, max_size)
                node_list.update(children)
            
        return node_list

    print('finding terminal nodes for {}.'.format(args.maxsize))
    nodes = find_terminal_nodes(top_node_id, args.maxsize)

    nodes = { n for n in nodes if node_sizes[n] > args.minsize }
    
    print('{} nodes between {} and {} in k-mer size'.format(len(nodes), args.minsize, args.maxsize))
    print('containing {} level1 nodes of {} total'.format(len(find_shadow(nodes, dag)), len(layer1_to_cdbg)))

    node_kmers = sum([ node_sizes[n] for n in nodes ])
    print('containing {} kmers of {} total'.format(node_kmers, node_sizes[top_node_id]))

    ### now build cdbg -> group ID

    cdbg_to_group = {}
    for n in nodes:
        shadow = find_shadow([n], dag)
        for level1_node in shadow:
            for cdbg_id in layer1_to_cdbg[level1_node]:
                assert cdbg_id not in cdbg_to_group
                cdbg_to_group[cdbg_id] = n

    ### record group info
    group_info = {}
    for n in nodes:
        group_info[n] = sourmash_lib.MinHash(n=0, ksize=4, scaled=1, track_abundance=1)

    ### aaaaaand iterate over contigs.
    for record_n, record in enumerate(screed.open(contigs)):
        if record_n % 10000 == 0:
            print('...', record_n, end='\r')
        cdbg_id = int(record.name)
        group_id = cdbg_to_group.get(cdbg_id)
        
        if group_id is not None:
            # keep/measure!
            mh = group_info[group_id]
            mh.add_sequence(record.sequence, True)

    D = numpy.zeros((len(group_info), len(group_info)))
    for i, (group_id, mh1) in enumerate(group_info.items()):
        for j, (group_id, mh2) in enumerate(group_info.items()):
            if i <= j:
                D[i][j] = mh1.similarity(mh2)
                D[j][i] = D[i][j]

    with open(args.output, 'wb') as fp:
        numpy.save(fp, D)

    with open(args.output + '.labels.txt', 'wt') as fp:
        for group_id in group_info:
            fp.write('{}\n'.format(group_id))


if __name__ == '__main__':
    main()
