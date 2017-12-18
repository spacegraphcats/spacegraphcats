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
import collections

import screed

import sourmash_lib
from sourmash_lib import MinHash

from .extract_reads import collect_shadow, collect_frontier
from .search_utils import (get_minhashdb_name, build_queries_for_seeds)
from spacegraphcats.logging import log
from . import search_utils
from .search_utils import (load_dag, load_layer1_to_cdbg)
from .frontier_search import find_shadow
from .extract_contigs_conn import assign_color, get_nodes


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query_seqs', help='query sequences')
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('overhead', help='\% of overhead', type=float)
    p.add_argument('-o', '--output', type=argparse.FileType('wt'))
    p.add_argument('--no-empty', action='store_true')
    p.add_argument('--purgatory', action='store_true')
    p.add_argument('--fullstats', action='store_true')
    p.add_argument('-k', '--ksize', default=31, type=int,
                   help='k-mer size (default: 31)')
    p.add_argument('--scaled', default=1000, type=float)
    p.add_argument('--seeds', default="43", type=str)
    p.add_argument('-v', '--verbose', action='store_true')
    p.add_argument('--diffuse-radius', type=int, default=0)

    args = p.parse_args()

    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')
    cdbgfile = os.path.join(args.catlas_prefix, 'cdbg.gxt')

    ##
    with open(cdbgfile, 'rt') as infp:
        cdbg = collections.defaultdict(set)
        num_nodes = int(next(infp))

        for line in infp:
            u, v = line.split()
            u, v = int(u), int(v)
            cdbg[u].add(v)
            cdbg[v].add(u)

    top_node_id, dag, cdbg_to_catlas = search_utils.load_just_dag(catlas)
    print('loaded catlas DAG with {} nodes'.format(len(dag)))

    contigs = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

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

    # build a query sig for each seed.
    print('building query signatures for {} seeds.'.format(len(seeds)))
    seed_queries = build_queries_for_seeds(seeds, ksize, scaled,
                                           args.query_seqs)

    # gather results of all queries across all seeds into total_shadow
    # gather results of all queries across all seeds into cdbg_shadow
    total_frontier = collect_frontier(seed_queries, dag, top_node_id,
                                      minhash_db_list, verbose=args.verbose,
                                      overhead=args.overhead)


    # load mapping between dom nodes and cDBG/graph nodes:
    layer1_to_cdbg = load_layer1_to_cdbg(cdbg_to_catlas, domfile)
    print('loaded {} layer 1 catlas nodes'.format(len(layer1_to_cdbg)))

    # assign colors
    total_cdbg_shadow = set()
    for x in total_frontier:
        cdbg_shadow = set()
        for y in find_shadow([x], dag):
            cdbg_shadow.update(layer1_to_cdbg.get(y, set()))
        assign_color(cdbg_shadow)
        total_cdbg_shadow.update(cdbg_shadow)

    next_layer = set()
    this_layer = total_cdbg_shadow

    #print(-1, group_max, color_max, len(set(group_to_color.values())))
    for dr in range(args.diffuse_radius):
        for u in this_layer:
            for v in cdbg.get(u, ()):
                next_layer.add(v)
                assign_color([u, v])

        this_layer = next_layer
        next_layer = set()

        #print(dr, group_max, color_max, len(set(group_to_color.values())))

    cdbg_shadow = set(get_nodes())

    # done with main loop! now extract contigs corresponding to union of
    # query results.
    print('------')
    print('done searching! now -> extracting contigs.')

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

    # track minhash of retrieved reads using original query minhash:
    reads_minhash = query_sig.minhash.copy_and_clear()

    print('loading contigs...')
    n = 0
    for record in screed.open(contigs):
        if n % 10000 == 0:
            offset_f = n / len(cdbg_shadow)
            print('...at n {} ({:.1f}% of shadow)'.format(n, offset_f * 100),
                  end='\r')

        contig_id = int(record.name)
        if contig_id not in cdbg_shadow:
            continue

        n += 1

        # track retrieved minhash
        reads_minhash.add_sequence(str(record.sequence), True)

        # output if desired!
        if args.output:
            args.output.write('>{}\n{}\n'.format(record.name, record.sequence))

        total_bp += len(record.sequence)
        total_seqs += 1

    print('')
    print('fetched {} contigs, {} bp matching combined frontiers.'.format(total_seqs, total_bp))

    print('orthogonal inclusion by retrieved contigs:', query_sig.minhash.contained_by(reads_minhash))

    end = time.time()
    print('total contig retrieval time (including database query): {:.2f}s'.format(end - start))

    sys.exit(0)


if __name__ == '__main__':
    main()
