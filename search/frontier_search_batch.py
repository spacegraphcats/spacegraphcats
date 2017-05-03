#! /usr/bin/env python
import os
import sys
import argparse
from sourmash_lib import MinHash
from collections import defaultdict
from spacegraphcats.catlas import CAtlas
import time
import sourmash_lib
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import search_minhashes, SigLeaf
from sourmash_lib import signature
import screed
import pickle
from typing import List

from .search_catlas_with_minhash import load_dag, load_minhash
from .frontier_search import frontier_search


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
    minhash_dir = os.path.join(args.catlas_prefix, basename + '.minhashes')

    for filename in args.query_sigs:
        query_sig = sourmash_lib.signature.load_signatures(filename)
        query_sig = list(query_sig)[0]
        print('loaded query sig {}'.format(query_sig.name()), file=sys.stderr)
        
        frontier, num_leaves, frontier_mh = frontier_search(query_sig, top_node_id, dag, minhash_dir, args.overhead)

        containment = query_mh.containment(frontier_mh)
        similarity = query_mh.similarity(frontier_mh)

        print(filename, containment, similarity)
    
    
if __name__ == '__main__':
    main()
