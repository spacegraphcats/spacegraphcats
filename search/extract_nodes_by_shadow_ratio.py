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
import pandas

import screed
from . import search_utils
from .frontier_search import find_shadow


def main(args=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('output')
    p.add_argument('--minsize', type=float, default=100)
    p.add_argument('--maxsize', type=float, default=10000)
    p.add_argument('-k', '--ksize', default=31, type=int,
                        help='k-mer size (default: 31)')
    args = p.parse_args(args)

    print('minsize: {:g}'.format(args.minsize))
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

    # calculate the cDBG shadow sizes for each catlas node.
    x = []
    for (node_id, level) in dag_levels.items():
        x.append((level, node_id))
    x.sort()

    node_shadow_sizes = {}
    for level, node_id in x:
        if level == 1:
            node_shadow_sizes[node_id] = len(layer1_to_cdbg[node_id])
        else:
            sub_size = 0
            for child_id in dag[node_id]:
                sub_size += node_shadow_sizes[child_id]
            node_shadow_sizes[node_id] = sub_size

    # ...and load cdbg node sizes
    print('loading contig size info')
    cdbg_kmer_sizes = search_utils.load_cdbg_size_info(args.catlas_prefix)

    # decorate catlas with cdbg node sizes underneath them
    x = []
    for (node_id, level) in dag_levels.items():
        x.append((level, node_id))
    x.sort()

    node_kmer_sizes = {}
    for level, node_id in x:
        if level == 1:
            total_kmers = 0
            for cdbg_node in layer1_to_cdbg.get(node_id):
                total_kmers += cdbg_kmer_sizes[cdbg_node]

            node_kmer_sizes[node_id] = total_kmers
        else:
            sub_size = 0
            for child_id in dag[node_id]:
                sub_size += node_kmer_sizes[child_id]
            node_kmer_sizes[node_id] = sub_size


    ### ok, the real work: look at articulation of cDBG graph.

    # find highest nodes with kmer size less than given max_size
    def find_terminal_nodes(node_id, max_size):
        node_list = set()
        for sub_id in dag[node_id]:
            # shadow size
            size = node_kmer_sizes[sub_id]

            if size < max_size:
                node_list.add(sub_id)
            else:
                children = find_terminal_nodes(sub_id, max_size)
                node_list.update(children)

        return node_list

    print('finding terminal nodes for {}.'.format(args.maxsize))

    terminal = find_terminal_nodes(top_node_id, args.maxsize)
    print('...got {}'.format(len(terminal)))
    terminal = { n for n in terminal if node_kmer_sizes[n] > args.minsize }
    print('...down to {} between {} and {} in size.'.format(len(terminal),
                                                            args.minsize,
                                                            args.maxsize))

    # now, go through and calculate total k-mers underneath.
    x = []
    for node_id in terminal:
        # calculate: how many k-mers per cDBG node?
        kmer_size = node_kmer_sizes[node_id]
        shadow_size = node_shadow_sizes[node_id]

        ratio = math.log(kmer_size, 2) - math.log(shadow_size, 2)

        # track basic info
        x.append((ratio, node_id, shadow_size, kmer_size))

    print('terminal node stats for maxsize: {:g}'.format(args.maxsize))
    print('n tnodes:', len(terminal))
    print('total k-mers:', node_kmer_sizes[top_node_id])

    x.sort(reverse=True)
    for (k, v, a, b) in x:
        print('ratio: {:.3f}'.format(2**k), '/ shadow size:', a, '/ kmers:', b)

    if 0:
        # keep only the last 10 (just for diagnostics/evaluation).
        keep_num = 10
        print('keeping last {} for examination.'.format(keep_num))
        keep_terminal = { n for (k, n, a, b) in x[-keep_num:] }
    elif 1:
        # keep the last 500kb for examination
        keep_sum_kmer = 500000
        sofar = 0
        keep_terminal = set()
        for (k, v, a, b) in reversed(x):
            sofar += b
            if sofar > keep_sum_kmer:
                break
            keep_terminal.add(v)

        print('keeping last {} k-mers worth of nodes for examination.'.format(sofar))

    # build cDBG shadow ID list.
    cdbg_shadow = set()
    terminal_shadow = find_shadow(keep_terminal, dag)
    for x in terminal_shadow:
        cdbg_shadow.update(layer1_to_cdbg.get(x))

    #### extract contigs
    print('extracting contigs & building a sourmash signature')
    contigs = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

    # track results as signature
    contigs_mh = sourmash_lib.MinHash(n=0, ksize=args.ksize, scaled=1000)

    total_bp = 0
    total_seqs = 0

    outfp = open(args.output, 'wt')
    for n, record in enumerate(screed.open(contigs)):
        if n and n % 10000 == 0:
            offset_f = total_seqs / len(cdbg_shadow)
            print('...at n {} ({:.1f}% of shadow)'.format(total_seqs,
                  offset_f * 100),
                  end='\r')

        # contig names == cDBG IDs
        contig_id = int(record.name)
        if contig_id not in cdbg_shadow:
            continue

        outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
        contigs_mh.add_sequence(record.sequence)

        # track retrieved sequences in a minhash
        total_bp += len(record.sequence)
        total_seqs += 1

    # done - got all contigs!
    print('')
    print('fetched {} contigs, {} bp.'.format(total_seqs, total_bp))

    print('wrote contigs to {}'.format(args.output))
    with open(args.output + '.sig', 'wt') as fp:
        ss = sourmash_lib.SourmashSignature(contigs_mh)
        sourmash_lib.save_signatures([ss], fp)


if __name__ == '__main__':
    main()
