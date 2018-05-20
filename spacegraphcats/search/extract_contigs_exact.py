#! /usr/bin/env python
"""
Retrieve the contigs for a list of cDBG nodes using exact/MPHF code.
"""
import argparse
import os
import sys
import gzip
import khmer
import collections
import bbhash
import numpy
from collections import defaultdict

import screed

from spacegraphcats.logging import log
from . import search_utils


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('query')
    p.add_argument('-k', '--ksize', default=31, type=int,
                   help='k-mer size (default: 31)')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'))
    p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args()

    contigs = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

    assert args.output, 'must specify -o'
    outfp = args.output
    outname = args.output.name

    # load k-mer MPHF index
    kmer_idx = search_utils.load_kmer_index(args.catlas_prefix)

    # build hashes for all the query k-mers
    print('loading query kmers...')
    bf = khmer.Nodetable(args.ksize, 1, 1)

    x = set()
    n = 0

    query_kmers = set()
    for record in screed.open(args.query):
        query_kmers.update(bf.get_kmer_hashes(record.sequence))

    cdbg_match_counts = kmer_idx.get_match_counts(query_kmers)

    f_found = sum(cdbg_match_counts.values()) / len(query_kmers)
    print('...done loading & counting query k-mers in cDBG.')
    print('containment: {:.1f}%'.format(f_found * 100))

    print('loading catlas...', end=' ')
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')
    top_node_id, dag, dag_up, dag_levels, catlas_to_cdbg = search_utils.load_dag(catlas)
    layer1_to_cdbg = search_utils.load_layer1_to_cdbg(catlas_to_cdbg, domfile)
    print('done.')

    cdbg_shadow = set()
    catlas_nodes = set()
    for layer1_node, cdbg_list in layer1_to_cdbg.items():
        for cdbg_node in cdbg_list:
            if cdbg_match_counts.get(cdbg_node, 0):
                # keep catlas node & all associated cdbg nodes
                catlas_nodes.add(layer1_node)
                cdbg_shadow.update(cdbg_list)
                break
    
    total_bp = 0
    total_seqs = 0

    # do some nodelist post-processing & output a response curve
    response_curve_filename = os.path.basename(outname) + '.response.txt'
    search_utils.output_response_curve(response_curve_filename,
                                       cdbg_match_counts, kmer_idx, layer1_to_cdbg)

    print('extracting contigs to {}.'.format(outname))
    for n, record in enumerate(screed.open(contigs)):
        if n % 10000 == 0:
            offset_f = total_seqs / len(cdbg_shadow)
            print('...at n {} ({:.1f}% of shadow)'.format(total_seqs,
                                                          offset_f * 100),
                  end='\r')

        contig_id = int(record.name)
        if contig_id not in cdbg_shadow:
            continue

        outfp.write('>{}\n{}\n'.format(record.name, record.sequence))

        total_bp += len(record.sequence)
        total_seqs += 1

    print('')
    print('fetched {} contigs, {} bp matching node list.'.format(total_seqs, total_bp))

    sys.exit(0)


if __name__ == '__main__':
    main()
