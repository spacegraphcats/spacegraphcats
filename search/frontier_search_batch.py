#! /usr/bin/env python
import argparse
import os
import pickle
import sys
import time
from collections import defaultdict
import csv
import leveldb

import screed
import sourmash_lib
from sourmash_lib import MinHash, signature
from sourmash_lib.sourmash_args import load_query_signature
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import SigLeaf, search_minhashes
from typing import List

from .frontier_search import frontier_search, compute_overhead
from .search_catlas_with_minhash import load_dag, load_minhash


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('query_sigs', help='query minhash list', nargs='+')
    p.add_argument('--overhead', help='\% of overhead', type=float,
                   default=0.1)
    p.add_argument('--purgatory', action='store_true')
    p.add_argument('-o', '--output', type=argparse.FileType('w'))
    p.add_argument('-k', '--ksize', default=None, type=int,
                        help='k-mer size (default: 31)')
    args = p.parse_args()

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')

    # load catlas DAG
    top_node_id, dag, dag_levels = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))
    db_path = os.path.join(args.catlas_prefix, 'minhashes.db')

    w = None
    if args.output:
        w = csv.writer(args.output)

    db = leveldb.LevelDB(os.path.join(args.catlas_prefix, 'minhashes.db'))

    for filename in args.query_sigs:
        query_sig = load_query_signature(filename, select_ksize=args.ksize,
                                         select_moltype='DNA')
        print('loaded query sig {}'.format(query_sig.name()), file=sys.stderr)
        
        frontier, _, _, frontier_mh = frontier_search(query_sig, top_node_id, dag, db, args.overhead, True, args.purgatory)

        query_mh = query_sig.minhash
        query_mh = query_mh.downsample_max_hash(frontier_mh)
        frontier_mh = frontier_mh.downsample_max_hash(query_mh)

        containment = query_mh.contained_by(frontier_mh)
        similarity = query_mh.similarity(frontier_mh)
        overhead = compute_overhead(frontier_mh, query_mh)

        print(filename, containment, similarity)
        if args.output:
            w.writerow([filename,containment,similarity,overhead])
    
    
if __name__ == '__main__':
    main()
