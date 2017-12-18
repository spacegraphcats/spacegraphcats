#! /usr/bin/env python
"""
Do a frontier search and retrieve contigs that match to the cDBG nodes
in the frontier.
"""
import argparse
import sys
import os
import sys
import gc
import csv
import traceback
import gzip

import leveldb
import screed

import sourmash_lib
from sourmash_lib import MinHash
from sourmash_lib.sourmash_args import load_query_signature

from .extract_reads import collect_frontier
from .search_utils import (get_minhashdb_name, get_reads_by_cdbg,
                           build_queries_for_seeds, parse_seeds_arg)
from spacegraphcats.logging import log
from search.frontier_search import (frontier_search, compute_overhead, find_shadow)
from . import search_utils
from .search_utils import (load_dag, load_layer1_to_cdbg, load_minhash)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('--overhead', help='\% of overhead', type=float, default=0.0)
    p.add_argument('output')
    p.add_argument('--query', help='query sequences', nargs='+')
    p.add_argument('--no-empty', action='store_true')
    p.add_argument('--purgatory', action='store_true')
    p.add_argument('--fullstats', action='store_true')
    p.add_argument('-k', '--ksize', default=31, type=int,
                   help='k-mer size (default: 31)')
    p.add_argument('--scaled', default=1000, type=float)
    p.add_argument('--seeds', default="43", type=str)

    args = p.parse_args()

    for filename in args.query:
        if not os.path.exists(filename):
            print('query seq file {} does not exist.'.format(filename))
            sys.exit(-1)

    try:
        os.mkdir(args.output)
    except OSError:
        pass
    if not os.path.isdir(args.output):
        print('output {} is not a directory'.format(args.output))
        sys.exit(-1)

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # load catlas DAG
    top_node_id, dag, dag_up, dag_levels, catlas_to_cdbg = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    # load mapping between dom nodes and cDBG/graph nodes:
    layer1_to_cdbg = load_layer1_to_cdbg(catlas_to_cdbg, domfile)
    print('loaded {} layer 1 catlas nodes'.format(len(layer1_to_cdbg)))

    # load contigs DB
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

    # record command line
    with open(os.path.join(args.output, 'command.txt'), 'wt') as fp:
        fp.write(str(sys.argv))

    # output results here:
    csvoutfp = open(os.path.join(args.output, 'results.csv'), 'wt')
    csv_writer = csv.writer(csvoutfp)
    csv_writer.writerow(['query', 'containment', 'similarity', 'bp', 'reads', 'n_seeds', 'ksize', 'scaled'])

    for query in args.query:
        try:
            print('------')
            print('QUERY FILE:', query)

            # build a query sig for each seed.
            seed_queries = build_queries_for_seeds(seeds, ksize, scaled, query)

            # gather results of all queries across all seeds
            total_frontier = collect_frontier(seed_queries, dag, top_node_id,
                                              minhash_db_list,
                                              overhead=args.overhead)

            total_shadow = find_shadow(total_frontier, dag)

            cdbg_shadow = set()
            for x in total_shadow:
                cdbg_shadow.update(layer1_to_cdbg.get(x))

            # done with main loop! now extract contigs corresponding to union
            # of query results.
            print('done searching! now -> extracting contigs.')

            # build check MinHash
            query_mh = MinHash(0, ksize, scaled=scaled, seed=42)

            name = None
            for record in screed.open(query):
                if not name:
                    name = record.name
                query_mh.add_sequence(record.sequence, True)
            query_sig = sourmash_lib.SourmashSignature(query_mh, name=name,
                                                       filename=query)

            total_bp = 0
            total_seqs = 0
            output_seqs = 0

            # track minhash of retrieved contigs using original query minhash:
            contigs_minhash = query_sig.minhash.copy_and_clear()

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
                contigs_minhash.add_sequence(str(record.sequence), True)

                total_bp += len(record.sequence)
                total_seqs += 1


            print('')
            print('fetched {} contigs, {} bp matching combined frontiers.'.format(total_seqs, total_bp))

            containment = query_sig.minhash.contained_by(contigs_minhash)
            similarity = query_sig.minhash.similarity(contigs_minhash)
            print('query inclusion by retrieved contigs:', containment)
            print('query similarity to retrieved contigs:', similarity)

            num_seeds = len(seeds)

            csv_writer.writerow([query, containment, similarity, total_bp, total_seqs, num_seeds, ksize, scaled])
            csvoutfp.flush()

            # write out signature from retrieved contigs.
            sig_filename = os.path.basename(query) + '.contigs.sig'
            with open(os.path.join(args.output, sig_filename), 'wt') as fp:
                ss = sourmash_lib.SourmashSignature(contigs_minhash,
                                                    name=sig_filename)
                sourmash_lib.save_signatures([ss], fp)

            # write out cDBG IDs
            cdbg_listname = os.path.basename(query) + '.cdbg_ids.txt.gz'
            with gzip.open(os.path.join(args.output, cdbg_listname), 'wt') as fp:
                fp.write("\n".join([ str(x) for x in cdbg_shadow ]))

            # write out frontier nodes by seed
            frontier_listname = os.path.basename(query) + '.frontier.txt.gz'
            with gzip.open(os.path.join(args.output, frontier_listname), 'wt') as fp:
                for node, seedlist in total_frontier.items():
                    fp.write('{},{}\n'.format(node,
                                              " ".join([ str(x) for x in seedlist ])))

        except:
            traceback.print_exc()

    ## end main loop!

    sys.exit(0)


if __name__ == '__main__':
    main()
