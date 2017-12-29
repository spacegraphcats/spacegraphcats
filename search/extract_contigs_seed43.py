#! /usr/bin/env python
"""
Retrieve the contigs for a list of cDBG nodes.  Consumes the output of
extract_nodes_by_query to get the list of nodes.
"""
import argparse
import os
import sys
import gzip
import khmer
import sourmash_lib

import screed

from spacegraphcats.logging import log
from . import search_utils
from .frontier_search import load_minhash


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('query')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'))
    p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args()

    contigs = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

    assert args.output, 'must specify -o'
    outfp = args.output
    outname = args.output.name

    query_mh = sourmash_lib.MinHash(n=0, ksize=31, scaled=1000, seed=43)
    for record in screed.open(args.query):
        query_mh.add_sequence(record.sequence)

    query_mins = query_mh.get_mins()

    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')
    top_node_id, dag, dag_up, dag_levels, catlas_to_cdbg = search_utils.load_dag(catlas)
    layer1_to_cdbg = search_utils.load_layer1_to_cdbg(catlas_to_cdbg, domfile)

    minhashdbfile = os.path.join(args.catlas_prefix, 'minhashes.db.k31.s1000.abund0.seed43')
    minhashdb = search_utils.MinhashSqlDB(minhashdbfile)

    n_empty = 0
    n_taken = 0
    n_iter = 0
    cdbg_shadow = set()
    for layer1_node, cdbg_list in layer1_to_cdbg.items():
        n_iter += 1
        assert dag_levels.get(layer1_node) == 1
        banded_mh = load_minhash(layer1_node, minhashdb)
        if banded_mh:
            banded_mins = set(banded_mh.get_mins())
            if banded_mins.intersection(query_mins):
                cdbg_shadow.update(cdbg_list)
                n_taken += 1
        else:
            n_empty += 1

    total_bp = 0
    total_seqs = 0

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

    print('fetched {} contigs, {} bp matching node list.'.format(total_seqs, total_bp))

    sys.exit(0)


if __name__ == '__main__':
    main()
