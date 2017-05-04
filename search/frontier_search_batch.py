#! /usr/bin/env python
import argparse
import os
import pickle
import sys
import time
from collections import defaultdict

import screed
import sourmash_lib
from sourmash_lib import MinHash, signature
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import SigLeaf, search_minhashes
from typing import List

from .frontier_search import frontier_search
from .search_catlas_with_minhash import load_dag, load_minhash


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('query_sigs', help='query minhash list', nargs='+')
    p.add_argument('--overhead', help='\% of overhead', type=float,
                   default=0.1)
    args = p.parse_args()

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, basename + '.catlas')

    # load catlas DAG
    top_node_id, dag, dag_levels = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))
    db_path = os.path.basename(args.catlas_prefix) + '/{}.db'.format(basename)

    for filename in args.query_sigs:
        query_sig = sourmash_lib.signature.load_signatures(filename)
        query_sig = list(query_sig)[0]
        print('loaded query sig {}'.format(query_sig.name()), file=sys.stderr)
        
        frontier, _, _, frontier_mh = frontier_search(query_sig, top_node_id, dag, db_path, args.overhead)

        containment = query_mh.containment(frontier_mh)
        similarity = query_mh.similarity(frontier_mh)

        print(filename, containment, similarity)
    
    
if __name__ == '__main__':
    main()
