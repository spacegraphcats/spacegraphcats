#! /usr/bin/env python
import argparse
import os
import pickle
import sys
import time
from collections import defaultdict
import csv
import leveldb
import gzip
import sqlite3

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
from . import bgzf
from . import search_utils
from .search_utils import get_minhashdb_name


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('readsfile')
    p.add_argument('labeled_reads_sqlite')
    p.add_argument('query_sigs', help='query minhash list', nargs='+')
    p.add_argument('--overhead', help='\% of overhead', type=float,
                   default=0.1)
    p.add_argument('--purgatory', action='store_true')
    p.add_argument('-o', '--output', type=argparse.FileType('w'))
    p.add_argument('-k', '--ksize', default=31, type=int,
                        help='k-mer size (default: 31)')
    p.add_argument('--savedir', default=None)
    args = p.parse_args()

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # load catlas DAG
    top_node_id, dag, dag_up, dag_levels = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

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

    db_path = get_minhashdb_name(args.catlas_prefix, args.ksize, 0, 0)
    if not db_path:
        print('** ERROR, minhash DB does not exist for {}'.format(args.ksize),
              file=sys.stderr)
        sys.exit(-1)
    minhash_db = leveldb.LevelDB(db_path)

    # sql
    dbfilename = args.labeled_reads_sqlite
    assert os.path.exists(dbfilename), 'sqlite file {} does not exist'.format(dbfilename)
    readsdb = sqlite3.connect(dbfilename)

    for filename in args.query_sigs:
        query_sig = load_query_signature(filename, args.ksize,
                                         select_moltype='DNA')
        print('loaded query sig {}'.format(query_sig.name()), file=sys.stderr)

        try:
            frontier, _, _, frontier_mh = frontier_search(query_sig, top_node_id, dag, minhash_db, args.overhead, True, args.purgatory)
        except:
            traceback.print_exc()
            continue
        
        print("removing empty frontier nodes...")
        nonempty_frontier = []

        for node in frontier:
            mh = load_minhash(node, minhash_db)
            if mh and len(mh.get_mins()) > 0:
                nonempty_frontier.append(node)
        print("...went from {} to {}".format(len(frontier), len(nonempty_frontier)))
        frontier = nonempty_frontier

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
            args.output.flush()

        if args.savedir:
            assert args.output
            outsig = os.path.basename(filename) + '.frontier'
            outsig = os.path.join(args.savedir, outsig)

            ## save frontier signature ##

            with open(outsig, 'w') as fp:
                sig = signature.SourmashSignature(frontier_mh,
                                                  name='frontier o={:1.2f} {}'.format(args.overhead, os.path.basename(filename)))
                sourmash_lib.signature.save_signatures([sig], fp)

            ## save reads ##

            cursor = readsdb.cursor()

            total_bp = 0
            watermark_size = 1e7
            watermark = watermark_size
            total_seqs = 0
            output_seqs = 0

            outreads = os.path.basename(filename) + '.reads.fa.gz'
            outreads = os.path.join(args.savedir, outreads)

            outfp = gzip.open(outreads, 'wt')
            reads_grabber = search_utils.GrabBGZF_Random(args.readsfile)

            ## get last offset:
            last_offset = search_utils.sqlite_get_max_offset(cursor)

            print('running sqlite query...')
            for n, offset in enumerate(search_utils.sqlite_get_offsets(cursor, cdbg_shadow)):

                if n % 10000 == 0:
                    print('...at n {} ({:.1f}% of file) - {} seqs'.format(n, offset /  last_offset * 100, total_seqs), end='\r')
                record, xx = reads_grabber.get_sequence_at(offset)
                assert xx == offset, (xx, offset)

                name = record.name
                sequence = record.sequence

                total_bp += len(sequence)
                total_seqs += 1

                outfp.write('>{}\n{}\n'.format(name, sequence))

            print('')
            print('fetched {} reads, {} bp matching frontier.'.format(total_seqs, total_bp))
            outfp.close()
                
    
if __name__ == '__main__':
    main()
