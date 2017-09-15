#! /usr/bin/env python
import argparse
import os
import sys
import leveldb
from copy import copy
import collections

import sourmash_lib
from sourmash_lib import MinHash, signature
from sourmash_lib.sourmash_args import load_query_signature
from typing import Dict, List, Set, Union, Tuple
from pickle import load
import screed
import sqlite3

from .memoize import memoize
from .search_catlas_with_minhash import load_dag, load_minhash
from spacegraphcats.logging import log
from search.frontier_search import (frontier_search, compute_overhead, find_shadow)
from .search_utils import (load_dag, load_layer0_to_cdbg)
import khmer.utils
from . import bgzf


def read_bgzf(reader):
    from screed.openscreed import fastq_iter, fasta_iter
    ch = reader.read(1)
    if ch == '>':
        iter_fn = fasta_iter
    elif ch == '@':
        iter_fn = fastq_iter
    else:
        raise Exception('unknown start chr {}'.format(ch))

    reader.seek(0)

    last_pos = reader.tell()
    for record in iter_fn(reader):
        yield record, last_pos
        last_pos = reader.tell()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query_sig', help='query minhash')
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('overhead', help='\% of overhead', type=float)
    p.add_argument('readsfile')
    p.add_argument('labeled_reads_sqlite')
    p.add_argument('output')
    p.add_argument('--no-empty', action='store_true')
    p.add_argument('--purgatory', action='store_true')
    p.add_argument('--fullstats', action='store_true')
    p.add_argument('-k', '--ksize', default=None, type=int,
                        help='k-mer size (default: 31)')

    args = p.parse_args()

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # load catlas DAG
    top_node_id, dag, dag_up, dag_levels = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    db_path = os.path.join(args.catlas_prefix, 'minhashes.db')
    minhash_db = leveldb.LevelDB(db_path)

    # load mapping between dom nodes and cDBG/graph nodes:
    layer0_to_cdbg = load_layer0_to_cdbg(catlas, domfile)
    print('loaded {} layer 0 catlas nodes'.format(len(layer0_to_cdbg)))
    x = set()
    for v in layer0_to_cdbg.values():
        x.update(v)
    print('...corresponding to {} cDBG nodes.'.format(len(x)))

    # load query MinHash
    query_sig = load_query_signature(args.query_sig, select_ksize=args.ksize,
                                     select_moltype='DNA')
    print('loaded query sig {}'.format(query_sig.name()))

    frontier, num_leaves, num_empty, frontier_mh = frontier_search(query_sig, top_node_id, dag, minhash_db, args.overhead, not args.no_empty, args.purgatory)

    top_mh = load_minhash(top_node_id, minhash_db)
    query_mh = query_sig.minhash.downsample_max_hash(top_mh)
    top_mh = top_mh.downsample_max_hash(query_sig.minhash)
    print("Root containment: {}".format(query_mh.contained_by(top_mh)))
    print("Root similarity: {}".format(query_mh.similarity(top_mh)))

    print("Containment of frontier: {}".format(query_mh.contained_by(frontier_mh)))
    print("Similarity of frontier: {}".format(query_mh.similarity(frontier_mh)))
    print("Size of frontier: {} of {} ({:.3}%)".format(len(frontier), len(dag), 100 * len(frontier) / len(dag)))
    print("Overhead of frontier: {}".format(compute_overhead(frontier_mh, query_mh)))
    print("Number of leaves in the frontier: {}".format(num_leaves))
    print("Number of empty catlas nodes in the frontier: {}".format(num_empty))
    print("")

    shadow = find_shadow(frontier, dag)

    print("Size of the frontier shadow: {}".format(len(shadow)))
    if len(shadow) == len(layer0_to_cdbg):
        print('\n*** WARNING: shadow is the entire graph! ***\n')

    query_size = len(query_sig.minhash.get_mins())
    query_bp = query_size * query_sig.minhash.scaled
    print("Size of query minhash: {} (est {:2.1e} bp)".\
              format(query_size, query_bp))
    minhash_size = len(frontier_mh.get_mins())
    minhash_bp = minhash_size * frontier_mh.scaled
    print("Size of frontier minhash: {} (est {:2.1e} bp); ratio {:.2f}".\
              format(minhash_size, minhash_bp, minhash_bp / query_bp))

    log(args.catlas_prefix, sys.argv)

    #### extract reads

    cdbg_shadow = set()
    for x in shadow:
        cdbg_shadow.update(layer0_to_cdbg.get(x))

    print('loading graph & labels/foo...')
    
    # actually don't need whole graph, just tags @CTB
    ng = khmer.Nodegraph(args.ksize, 1, 1)

    # sql
    dbfilename = args.labeled_reads_sqlite + '.sqlite'
    assert os.path.exists(dbfilename), 'sqlite file {} does not exist'.format(dbfilename)
    db = sqlite3.connect(dbfilename)
    cursor = db.cursor()

    total_bp = 0
    watermark_size = 1e7
    watermark = watermark_size
    no_tags = 0
    no_tags_bp = 0
    total_seqs = 0
    output_seqs = 0

    outfp = open(args.output, 'wt')
    reader = bgzf.BgzfReader(args.readsfile, 'rt')
    reads_iter = read_bgzf(reader)
    next(reads_iter)

    ## get last offset:
    cursor.execute('SELECT max(sequences.offset) FROM sequences')
    last_offset = list(cursor)[0][0]

    seen_seq = set()
    if 1:
        label = 0
        cursor.execute('CREATE TEMPORARY TABLE label_query (label_id INTEGER PRIMARY KEY);')
        for label in cdbg_shadow:
            cursor.execute('INSERT INTO label_query (label_id) VALUES (?)', (label,))

        print('running sqlite query...')

        cursor.execute('SELECT DISTINCT sequences.offset FROM sequences WHERE label in (SELECT label_id FROM label_query) ORDER BY offset')
        print('...fetching read offsets')

        # @CTB try sorting? => then can work with .gz file?
        for n, (offset,) in enumerate(cursor):
            if n % 10000 == 0:
                print('...at n {} ({:.1f}% of file)'.format(n, offset /  last_offset * 100), end='\r')
            reader.seek(offset)
            (record, xx) = next(reads_iter)

            name = record.name
            sequence = record.sequence

            assert xx == offset, (xx, offset)
            
            total_bp += len(sequence)
            total_seqs += 1

            outfp.write('>{}\n{}\n'.format(name, sequence))

        print('')
        print('fetched {} reads, {} bp matching frontier.'.format(total_seqs, total_bp))

    sys.exit(0)


if __name__ == '__main__':
    # import cProfile
    # cProfile.run('main()', 'search_stats')

    main()
