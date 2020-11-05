#! /usr/bin/env python
"""
Given a catlas, the output of index_cdbg_by_hashval, and a sourmash
LCA database, summarize taxonomic confusion of dom nodes based on the
taxonomic content of the cDBG nodes under them.

This is exploratory :)
"""
import sys
import argparse
import time
from collections import defaultdict, Counter
import csv
import pickle

import spacegraphcats
from spacegraphcats.utils.logging import notify, error, debug
from spacegraphcats.search import search_utils
from spacegraphcats.search import MPHF_KmerIndex, hash_sequence
from spacegraphcats.search.catlas import CAtlas
import screed
from sourmash import lca
from sourmash.lca import LCA_Database


def main():
    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('mh_index_picklefile', help='pickled hashval index')
    p.add_argument('lca_db')
    args = p.parse_args()

    # load catlas DAG
    catlas = CAtlas(args.catlas_prefix, load_sizefile=True)
    notify('loaded {} nodes from catlas {}', len(catlas), args.catlas_prefix)
    notify('loaded {} layer 1 catlas nodes', len(catlas.layer1_to_cdbg))

    ki_start = time.time()
    kmer_idx = MPHF_KmerIndex.from_catlas_directory(args.catlas_prefix)
    notify('loaded {} k-mers in index ({:.1f}s)',
           len(kmer_idx.mphf_to_kmer), time.time() - ki_start)

    # load picklefile
    with open(args.mh_index_picklefile, 'rb') as fp:
        hashval_to_contig_id = pickle.load(fp)
    notify('loaded {} hash value -> cdbg_id mappings from {}',
           len(hashval_to_contig_id), args.mh_index_picklefile)

    # load LCA DB
    db = LCA_Database.load(args.lca_db)

    ###

    # first, aggregate cDBG hashvals to dom node hashvals
    dom_hashvals = defaultdict(set)
    for hashval, cdbg_id in hashval_to_contig_id.items():
        dom_id = catlas.cdbg_to_layer1[cdbg_id]
        dom_hashvals[dom_id].add(hashval)

    # for each dom node, calculate LCA
    n_zero = 0
    n_one = 0
    n_many = 0
    dom_many_ids = set()
    for dom_id in catlas:
        if catlas.levels[dom_id] != 1:
            continue

        hashvals = dom_hashvals.get(dom_id)
        if not hashvals:
            n_zero += 1
        elif len(hashvals) == 1:
            n_one += 1
        else:
            n_many += 1
            dom_many_ids.add(dom_id)

    print(f'{n_zero} dom nodes have no sourmash hashes under them.')
    print(f'{n_one} dom nodes have exactly one sourmash hash under them.')
    print(f'{n_many} dom nodes have two or more sourmash hashes under them.')

    cnt = Counter()
    for dom_id in dom_many_ids:
        hashvals = dom_hashvals[dom_id]
        assert len(hashvals) > 1

        lins = set()
        for hashval in hashvals:
            lineages = db.get_lineage_assignments(hashval)
            lins.update(lineages)

        tree = lca.build_tree(lins)
        lca_lin, reason = lca.find_lca(tree)
        cnt[lca_lin[-1].rank] += 1

    print('')
    print(f'rank of dom node lca  count of dom nodes with that rank')
    print(f'--------------------  ---------------------------------')
    for rank, count in cnt.most_common():
        print(f'{rank}                  {count}')

    return 0


if __name__ == '__main__':
    sys.exit(main())
