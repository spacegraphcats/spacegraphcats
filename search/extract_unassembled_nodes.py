#! /usr/bin/env python
"""
Extract information about catlas nodes after splitting them into "assembled"
and "not".
"""
import argparse
import os
import sys
import collections
import math
import sourmash_lib
import khmer
import csv

import screed
from . import search_utils
from .frontier_search import find_shadow


def main(args=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('query')
    p.add_argument('output')
    p.add_argument('--threshold', default=0.0)
    p.add_argument('-k', '--ksize', default=31, type=int,
                   help='k-mer size (default: 31)')
    args = p.parse_args(args)

    print('threshold: {:.3f}'.format(args.threshold))

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
    print('decorating catlas with shadow size info.')
    node_shadow_sizes = search_utils.decorate_catlas_with_shadow_sizes(layer1_to_cdbg, dag, dag_levels)

    # ...and load cdbg node sizes
    print('loading contig size info')
    cdbg_kmer_sizes, cdbg_weighted_kmer_sizes = search_utils.load_cdbg_size_info(args.catlas_prefix)

    # decorate catlas with cdbg node sizes underneath them
    print('decorating catlas with contig size info.')
    node_kmer_sizes, node_weighted_kmer_sizes = search_utils.decorate_catlas_with_kmer_sizes(layer1_to_cdbg, dag, dag_levels, cdbg_kmer_sizes, cdbg_weighted_kmer_sizes)

    # load k-mer index, query, etc. etc.
    kmer_idx = search_utils.load_kmer_index(args.catlas_prefix)

    bf = khmer.Nodetable(args.ksize, 1, 1)

    query_kmers = set()
    for record in screed.open(args.query):
        query_kmers.update(bf.get_kmer_hashes(record.sequence))

    print('got {} k-mers from {}'.format(len(query_kmers), args.query))

    # construct dict cdbg_id -> # of query k-mers
    cdbg_match_counts = kmer_idx.get_match_counts(query_kmers)

    total_match_kmers = sum(cdbg_match_counts.values())
    f_found = total_match_kmers / len(query_kmers)
    print('=> containment: {:.1f}%'.format(f_found * 100))
    print('done loading & counting query k-mers in cDBG.')

    # calculate the cDBG matching k-mers sizes for each catlas node.
    catlas_match_counts = kmer_idx.build_catlas_match_counts(cdbg_match_counts, dag, dag_levels, layer1_to_cdbg)

    ### ok, the real work: find nodes that have low # of k-mers in the query.
    def find_unassembled_nodes(node_id, threshold=0.0):
        node_list = set()
        for sub_id in dag[node_id]:
            n_matched = catlas_match_counts.get(sub_id, 0)
            size = node_kmer_sizes[sub_id]

            f_assembled = n_matched / size

            if f_assembled <= threshold:
                node_list.add(sub_id)
            else:
                children = find_unassembled_nodes(sub_id, threshold)
                node_list.update(children)

        return node_list

    print('finding unassembled nodes for threshold {}.'.format(args.threshold))

    terminal = find_unassembled_nodes(top_node_id, args.threshold)
    sum_kmers = sum([ node_kmer_sizes[n] for n in terminal ])
    sum_match_kmers = sum([ catlas_match_counts.get(n, 0) for n in terminal ])
    print('...got {} nodes, representing {} k-mers ({} matched - {:.1f}'.format(len(terminal), sum_kmers, sum_match_kmers, sum_match_kmers / sum_kmers * 100))

    # now, go through all nodes and print out characteristics
    print('writing node info to {}'.format(args.output + '.csv'))
    with open(args.output + '.csv', 'wt') as fp:
        w = csv.writer(fp)

        w.writerow(['node_id', 'contained', 'n_kmers', 'n_weighted_kmers', 'shadow_size'])
        for n in terminal:
            f_contained = catlas_match_counts.get(n, 0) / node_kmer_sizes[n]
            w.writerow([str(n), str(f_contained), str(node_kmer_sizes[n]),
                        str(node_weighted_kmer_sizes[n]),
                        str(node_shadow_sizes[n])])

    # build cDBG shadow ID list.
    cdbg_shadow = set()
    terminal_shadow = find_shadow(terminal, dag)
    for x in terminal_shadow:
        cdbg_shadow.update(layer1_to_cdbg.get(x))

    #### extract contigs
    print('extracting contigs & building a sourmash signature')
    contigs = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

    # track results as signature
    contigs_mh = sourmash_lib.MinHash(n=0, ksize=args.ksize, scaled=1000)

    total_bp = 0
    total_seqs = 0

    print('writing contigs to {}'.format(args.output + '.fa'))
    outfp = open(args.output + '.fa', 'wt')
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

    print('writing sig to {}'.format(args.output + '.sig'))
    with open(args.output + '.sig', 'wt') as fp:
        ss = sourmash_lib.SourmashSignature(contigs_mh)
        sourmash_lib.save_signatures([ss], fp)



if __name__ == '__main__':
    main()
