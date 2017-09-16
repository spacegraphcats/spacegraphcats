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
import traceback

from .frontier_search import frontier_search, compute_overhead, find_shadow
from .search_catlas_with_minhash import load_dag, load_minhash
from .search_utils import (load_dag, load_layer0_to_cdbg)


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
    p.add_argument('--savedir', default=None)
    args = p.parse_args()

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # load catlas DAG
    top_node_id, dag, dag_up, dag_levels = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))
    db_path = os.path.join(args.catlas_prefix, 'minhashes.db')

    # load mapping between dom nodes and cDBG/graph nodes:
    layer0_to_cdbg = load_layer0_to_cdbg(catlas, domfile)
    print('loaded {} layer 0 catlas nodes'.format(len(layer0_to_cdbg)))
    x = set()
    for v in layer0_to_cdbg.values():
        x.update(v)
    print('...corresponding to {} cDBG nodes.'.format(len(x)))

    w = None
    if args.output:
        w = csv.writer(args.output)

    db = leveldb.LevelDB(os.path.join(args.catlas_prefix, 'minhashes.db'))

    for filename in args.query_sigs:
        query_sig = load_query_signature(filename, select_ksize=args.ksize,
                                         select_moltype='DNA')
        print('loaded query sig {}'.format(query_sig.name()), file=sys.stderr)

        try:
            frontier, _, _, frontier_mh = frontier_search(query_sig, top_node_id, dag, db, args.overhead, True, args.purgatory)
        except:
            traceback.print_exc()
            continue
        
        shadow = find_shadow(frontier, dag)
        cdbg_shadow = set()
        for x in shadow:
            cdbg_shadow.update(layer0_to_cdbg.get(x))

        query_mh = query_sig.minhash
        query_mh = query_mh.downsample_max_hash(frontier_mh)
        frontier_mh = frontier_mh.downsample_max_hash(query_mh)

        containment = query_mh.contained_by(frontier_mh)
        similarity = query_mh.similarity(frontier_mh)
        overhead = compute_overhead(frontier_mh, query_mh)

        print(filename, containment, similarity)
        if args.output:
            w.writerow([filename,containment,similarity,overhead])

        if args.savedir:
            assert args.output
            outsig = os.path.basename(filename) + '.frontier'
            outsig = os.path.join(args.savedir, outsig)

            with open(outsig, 'w') as fp:
                sig = signature.SourmashSignature('', frontier_mh,
                                                  name='frontier o={:1.2f} {}'.format(args.overhead, str(args.output)))
                sourmash_lib.signature.save_signatures([sig], fp)

            outshadow = os.path.basename(filename) + '.shadow'
            outshadow = os.path.join(args.savedir, outshadow)

            with open(outshadow, 'wb') as fp:
                pickle.dump(cdbg_shadow, fp)
                
    
if __name__ == '__main__':
    main()
