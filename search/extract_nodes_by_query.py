#! /usr/bin/env python
"""
Do a frontier search, and retrieve cDBG node IDs and MinHash signature for
the retrieved contigs.
"""
import argparse
import sys
import os
import sys
import gc
import csv
import traceback
import gzip
from collections import defaultdict
import time

import screed

import sourmash_lib
from sourmash_lib import MinHash
from sourmash_lib.sourmash_args import load_query_signature

from .search_utils import (get_minhashdb_name, get_reads_by_cdbg,
                           build_queries_for_seeds, parse_seeds_arg)
from .search_utils import MinhashSqlDB
from spacegraphcats.logging import log
from search.frontier_search import (frontier_search, compute_overhead, find_shadow)
from . import search_utils
from .search_utils import (load_dag, load_layer1_to_cdbg, load_minhash)




def collect_frontier(seed_queries, dag, top_node_id, minhash_db_list,
                     overhead=0.0, verbose=False):
    # gather results of all queries across all seeds into total_frontier
    total_frontier = defaultdict(set)

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

        for node in frontier:
            total_frontier[node].add(seed_query.minhash.seed)

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

    return total_frontier


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('--overhead', help='\% of overhead', type=float, default=0.0)
    p.add_argument('output')
    p.add_argument('--query', help='query sequences', nargs='+')
    p.add_argument('--no-empty', action='store_true')
    p.add_argument('-k', '--ksize', default=31, type=int,
                   help='k-mer size (default: 31)')
    p.add_argument('--scaled', default=1000, type=float)
    p.add_argument('--seeds', default="43", type=str)

    args = p.parse_args()

    # make sure all of the query sequences exist.
    for filename in args.query:
        if not os.path.exists(filename):
            print('query seq file {} does not exist.'.format(filename))
            sys.exit(-1)

    # create output directory if it doesn't exist.
    try:
        os.mkdir(args.output)
    except OSError:
        pass
    if not os.path.isdir(args.output):
        print('output {} is not a directory'.format(args.output))
        sys.exit(-1)

    # figure out catlas and domfile information.
    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # load catlas DAG
    top_node_id, dag, dag_up, dag_levels, catlas_to_cdbg = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    del dag_up
    del dag_levels

    gc.collect()

    # load mapping between dom nodes and cDBG/graph nodes:
    layer1_to_cdbg = load_layer1_to_cdbg(catlas_to_cdbg, domfile)
    print('loaded {} layer 1 catlas nodes'.format(len(layer1_to_cdbg)))

    del catlas_to_cdbg

    gc.collect()

    # find the contigs filename
    contigs = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

    # get a single ksize & scaled
    ksize = int(args.ksize)
    scaled = int(args.scaled)

    # we'll probably be using multiple seeds...
    seeds = search_utils.parse_seeds_arg(args.seeds)
    print('seeds: {} -> {}'.format(args.seeds, seeds))
    assert 42 not in seeds

    # make sure there's a minhash db for each seed 
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
        fp.write("\n")

    # output results.csv in the output directory:
    csvoutfp = open(os.path.join(args.output, 'results.csv'), 'wt')
    csv_writer = csv.writer(csvoutfp)
    csv_writer.writerow(['query', 'containment', 'similarity', 'bp', 'reads', 'n_seeds', 'ksize', 'scaled'])

    # iterate over each query, do the thing.
    for query in args.query:
        # ignore all the problems!
        try:
            print('------')
            print('QUERY FILE:', query)

            # build a query sig for each seed.
            seed_queries = build_queries_for_seeds(seeds, ksize, scaled, query)

            # gather results of all queries across all seeds
            total_frontier = collect_frontier(seed_queries, dag, top_node_id,
                                              minhash_db_list,
                                              overhead=args.overhead)

            # calculate level 1 nodes for this frontier in the catlas
            total_shadow = find_shadow(total_frontier, dag)

            # calculate associated cDBG nodes
            cdbg_shadow = set()
            for x in total_shadow:
                cdbg_shadow.update(layer1_to_cdbg.get(x))

            # done with main loop! now extract contigs using cDBG shadow
            # node list.
            print('done searching! now -> extracting contigs.')

            # track extracted info
            total_bp = 0
            total_seqs = 0

            # build check MinHash w/seed=42
            query_sig = build_queries_for_seeds([42], ksize, scaled, query)[0]

            # track minhash of retrieved contigs using original query minhash:
            contigs_minhash = query_sig.minhash.copy_and_clear()

            # walk through the contigs, retrieving.
            print('loading contigs...')
            for n, record in enumerate(screed.open(contigs)):
                if n and n % 10000 == 0:
                    offset_f = total_seqs / len(cdbg_shadow)
                    print('...at n {} ({:.1f}% of shadow)'.format(total_seqs,
                          offset_f * 100),
                          end='\r')

                # contig names == cDBG IDs
                contig_id = int(record.name)
                if contig_id not in cdbg_shadow:
                    continue

                # track retrieved sequences in a minhash
                contigs_minhash.add_sequence(str(record.sequence), True)

                total_bp += len(record.sequence)
                total_seqs += 1

            # done - got all contigs!
            print('')
            print('fetched {} contigs, {} bp matching combined frontiers.'.format(total_seqs, total_bp))

            # calculate summary values of extracted contigs
            containment = query_sig.minhash.contained_by(contigs_minhash)
            similarity = query_sig.minhash.similarity(contigs_minhash)
            print('query inclusion by retrieved contigs:', containment)
            print('query similarity to retrieved contigs:', similarity)

            num_seeds = len(seeds)

            # output to results.csv!
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

        except KeyboardInterrupt:
            raise
        except:
            traceback.print_exc()

    ## end main loop!

    sys.exit(0)


if __name__ == '__main__':
    main()
