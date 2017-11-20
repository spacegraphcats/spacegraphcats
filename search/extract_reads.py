#! /usr/bin/env python
"""
Do a frontier search and retrieve reads that match to the cDBG nodes
in the frontier.

Uses the output of label_cdbg to do its magic.

See also extract_contigs_by_frontier for a simpler version that just pulls
out linear segments.
"""
import argparse
import os
import sys
import time
import gc

import screed

import sourmash_lib
from sourmash_lib import MinHash
from sourmash_lib.sourmash_args import load_query_signature

from .search_utils import (get_minhashdb_name, get_reads_by_cdbg,
                           build_queries_for_seeds, MinhashSqlDB)
from spacegraphcats.logging import log
from search.frontier_search import (frontier_search, compute_overhead, find_shadow)
from . import search_utils
from .search_utils import (load_dag, load_layer1_to_cdbg, load_minhash)


def collect_shadow(seed_queries, dag, top_node_id, minhash_db_list,
                   overhead=0.0, verbose=False):
    # gather results of all queries across all seeds into cdbg_shadow
    total_shadow = set()

    # do queries!
    for seed_query, db_path in zip(seed_queries, minhash_db_list):
        start = time.time()
        print('loading minhashdb:', db_path)
        minhash_db = MinhashSqlDB(db_path)

        print('searching with seed={}'.format(seed_query.minhash.seed))
        frontier, num_leaves, num_empty, frontier_mh = \
          frontier_search(seed_query, top_node_id, dag, minhash_db,
                          overhead, False, False)

        top_mh = load_minhash(top_node_id, minhash_db)
        query_mh = seed_query.minhash.downsample_max_hash(top_mh)
        top_mh = top_mh.downsample_max_hash(seed_query.minhash)

        if verbose:
            print("Root containment: {}".format(query_mh.contained_by(top_mh)))
            print("Root similarity: {}".format(query_mh.similarity(top_mh)))

            print("Containment of frontier: {}".format(query_mh.contained_by(frontier_mh)))
            print("Similarity of frontier: {}".format(query_mh.similarity(frontier_mh)))
            print("Size of frontier: {} of {} ({:.3}%)".format(len(frontier), len(dag), 100 * len(frontier) / len(dag)))
            print("Overhead of frontier: {}".format(compute_overhead(frontier_mh, query_mh)))
            print("Number of leaves in the frontier: {}".format(num_leaves))
            print("Number of empty catlas nodes in the frontier: {}".format(num_empty))
            print("")

        if verbose:
            print("removing empty catlas nodes from the frontier...")
        nonempty_frontier = search_utils.remove_empty_catlas_nodes(frontier,
                                                                   minhash_db)

        if verbose:
            print("...went from {} to {}".format(len(frontier), len(nonempty_frontier)))
        frontier = nonempty_frontier

        shadow = find_shadow(frontier, dag)
        total_shadow.update(shadow)

        if verbose:
            print("Size of the frontier shadow: {}".format(len(shadow)))

        query_size = len(seed_query.minhash.get_mins())
        query_bp = query_size * seed_query.minhash.scaled
        if verbose:
            print("Size of query minhash: {} (est {:2.1e} bp)".\
                    format(query_size, query_bp))
        minhash_size = len(frontier_mh.get_mins())
        minhash_bp = minhash_size * frontier_mh.scaled
        if verbose:
            print("Size of frontier minhash: {} (est {:2.1e} bp); ratio {:.2f}".\
                  format(minhash_size, minhash_bp, minhash_bp / query_bp))

        end = time.time()
        print('query time: {:.1f}s'.format(end-start))

    return total_shadow


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query_seqs', help='query sequences')
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('overhead', help='\% of overhead', type=float, default=0.0)
    p.add_argument('readsfile')
    p.add_argument('labeled_reads_sqlite')
    p.add_argument('output')
    p.add_argument('--no-empty', action='store_true')
    p.add_argument('--purgatory', action='store_true')
    p.add_argument('--fullstats', action='store_true')
    p.add_argument('-k', '--ksize', default=31, type=int,
                   help='k-mer size (default: 31)')
    p.add_argument('--scaled', default=1000, type=float)
    p.add_argument('--seeds', default="43", type=str)
    p.add_argument('-v', '--verbose', action='store_true')

    args = p.parse_args()

    dbfilename = args.labeled_reads_sqlite
    if not os.path.exists(dbfilename):
        print('sqlite file {} does not exist'.format(dbfilename))
        sys.exit(-1)

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # load catlas DAG
    top_node_id, dag, cdbg_to_catlas = search_utils.load_just_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    # single ksize
    ksize = int(args.ksize)
    scaled = int(args.scaled)

    # multiple seeds
    seeds = search_utils.parse_seeds_arg(args.seeds)
    print('seeds: {} -> {}'.format(args.seeds, seeds))
    assert 42 not in seeds

    # locate a minhash db for each seed
    minhash_db_list = []
    for seed in seeds:
        db_path = get_minhashdb_name(args.catlas_prefix, ksize, scaled, 0,
                                     seed)
        if not db_path:
            print('** ERROR, minhash DB does not exist for k={} seed={}'.format(ksize, seed),
                  file=sys.stderr)
            sys.exit(-1)
        minhash_db_list.append(db_path)

    log(args.catlas_prefix, sys.argv)

    # build a query sig for each seed.
    print('building query signatures for {} seeds.'.format(len(seeds)))
    seed_queries = build_queries_for_seeds(seeds, ksize, scaled,
                                           args.query_seqs)

    # gather results of all queries across all seeds into cdbg_shadow
    total_shadow = collect_shadow(seed_queries, dag, top_node_id,
                                  minhash_db_list, verbose=args.verbose,
                                  overhead=args.overhead)

    del dag
    gc.collect()
    
    # load mapping between dom nodes and cDBG/graph nodes:
    layer1_to_cdbg = load_layer1_to_cdbg(cdbg_to_catlas, domfile)
    print('loaded {} layer 1 catlas nodes'.format(len(layer1_to_cdbg)))

    cdbg_shadow = set()
    for x in total_shadow:
        cdbg_shadow.update(layer1_to_cdbg.get(x))

    # done with main loop! now extract reads corresponding to union of
    # query results.
    print('------')
    print('done searching! now -> extracting reads.')

    # build check MinHash
    print('building check MinHash with orthogonal seed=42')
    query_mh = MinHash(0, ksize, scaled=scaled, seed=42)

    name = None
    for record in screed.open(args.query_seqs):
        if not name:
            name = record.name
        query_mh.add_sequence(record.sequence, True)
    print('loaded check sequences from {}'.format(args.query_seqs))
    query_sig = sourmash_lib.SourmashSignature(query_mh, name=name,
                                               filename=args.query_seqs)


    start = time.time()
    
    total_bp = 0
    total_seqs = 0
    output_seqs = 0

    # output sequences here:
    outfp = open(args.output, 'wt')

    # track minhash of retrieved reads using original query minhash:
    reads_minhash = query_sig.minhash.copy_and_clear()

    print('loading labels and querying for matching sequences..')
    reads_iter = get_reads_by_cdbg(dbfilename, args.readsfile, cdbg_shadow)
    for n, (record, offset_f) in enumerate(reads_iter):
        if n % 10000 == 0:
            print('...at n {} ({:.1f}% of file)'.format(n, offset_f * 100),
                  end='\r')

        # track retrieved minhash
        reads_minhash.add_sequence(record.sequence, True)

        # output!
        outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
        total_bp += len(record.sequence)
        total_seqs += 1

    print('')
    print('fetched {} reads, {} bp matching combined frontiers.'.format(total_seqs, total_bp))

    print('orthogonal inclusion by retrieved reads:', query_sig.minhash.contained_by(reads_minhash))

    end = time.time()
    print('total read retrieval time (including database query): {:.2f}s'.format(end - start))

    sys.exit(0)


if __name__ == '__main__':
    main()
