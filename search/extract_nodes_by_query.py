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
import khmer
import pickle

import screed

import sourmash_lib
from sourmash_lib import MinHash
from sourmash_lib.sourmash_args import load_query_signature
from sourmash_lib._minhash import hash_murmur

from .search_utils import get_reads_by_cdbg, load_kmer_index
from spacegraphcats.logging import log
from search.frontier_search import (frontier_search,
                                    frontier_search_exact,
                                    find_shadow,
                                    NoContainment)
from . import search_utils
from .search_utils import (load_dag, load_layer1_to_cdbg)


def build_query_mh_for_seed(seed, ksize, scaled, query_seq_file):
    mh = MinHash(0, ksize, scaled=scaled, seed=seed)

    name = None
    for record in screed.open(query_seq_file):
        if not name:
            name = record.name

        mh.add_sequence(record.sequence, True)

    return sourmash_lib.SourmashSignature(mh, name=name,
                                          filename=query_seq_file)


def collect_frontier(dag,
                     top_node_id,
                     node_sizes,
                     catlas_match_counts,
                     max_overhead,
                     min_containment,
                     verbose=False):
    start = time.time()

    # gather results into total_frontier
    total_frontier = defaultdict(set)

    # do queries!
    try:
        frontier, num_leaves, num_empty, frontier_mh = \
          frontier_search(top_node_id,
                          dag,
                          node_sizes,
                          catlas_match_counts,
                          max_overhead,
                          min_containment)
    except NoContainment:
        print('** WARNING: no containment!?')
        frontier = []

    # record which seed (always 0, here) contributed to which node
    for node in frontier:
        total_frontier[node].add(0)

    end = time.time()
    print('catlas query time: {:.1f}s'.format(end-start))

    return total_frontier


def collect_frontier_exact(dag,
                           top_node_id,
                           node_sizes,
                           catlas_match_counts,
                           overhead=0.0,
                           verbose=False):
    start = time.time()

    # gather results into total_frontier
    total_frontier = defaultdict(set)

    # do queries!
    try:
        frontier, num_leaves, num_empty, frontier_mh = \
          frontier_search_exact(top_node_id,
                                dag,
                                node_sizes,
                                catlas_match_counts,
                                overhead)

    except NoContainment:
        print('** WARNING: no containment!?')
        frontier = []
        num_leaves = 0
        num_empty = 0
        frontier_mh = seed_query.minhash.copy_and_clear()

    # record which seed (always 0, here) contributed to which node
    for node in frontier:
        total_frontier[node].add(0)

    end = time.time()
    print('catlas query time: {:.1f}s'.format(end-start))

    return total_frontier


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('--overhead', help='\% of overhead', type=float,
                   default=0.0)
    p.add_argument('output')
    p.add_argument('--min_containment', help="minimum containment",
                   type=float, default=1.0)
    p.add_argument('--max_overhead', help="largest overhead allowed",
                   type=float, default=1.0)
    p.add_argument('--query', help='query sequences', nargs='+')
    p.add_argument('--no-empty', action='store_true')
    p.add_argument('-k', '--ksize', default=31, type=int,
                   help='k-mer size (default: 31)')
    p.add_argument('--scaled', default=1000, type=float)
    p.add_argument('-v', '--verbose', action='store_true')

    args = p.parse_args(argv)

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

    # load mapping between dom nodes and cDBG/graph nodes:
    layer1_to_cdbg = load_layer1_to_cdbg(catlas_to_cdbg, domfile)
    print('loaded {} layer 1 catlas nodes'.format(len(layer1_to_cdbg)))

    # find the contigs filename
    contigs = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

    # ...and kmer index.
    ki_start = time.time()
    kmer_idx = load_kmer_index(args.catlas_prefix)
    print('loaded {} k-mers in index ({:.1f}s)'.format(
                                len(kmer_idx.mphf_to_kmer),
                                time.time() - ki_start))

    # calculate the k-mer sizes for each catlas node.
    node_sizes = kmer_idx.build_catlas_node_sizes(dag,
                                                  dag_levels,
                                                  layer1_to_cdbg)

    # get a single ksize & scaled
    ksize = int(args.ksize)
    scaled = int(args.scaled)

    # record command line
    with open(os.path.join(args.output, 'command.txt'), 'wt') as fp:
        fp.write(str(sys.argv))
        fp.write("\n")

    # output results.csv in the output directory:
    csvoutfp = open(os.path.join(args.output, 'results.csv'), 'wt')
    csv_writer = csv.writer(csvoutfp)
    csv_writer.writerow(['query', 'containment', 'similarity', 'bp', 'contigs',
                         'ksize', 'num_query_kmers', 'best_containment',
                         'cdbg_min_overhead', 'catlas_min_overhead'])

    # iterate over each query, do the thing.
    for query in args.query:
        # ignore all the problems!
        try:
            print('----')
            print('QUERY FILE:', query)
            start_time = time.time()

            # build hashes for all the query k-mers
            print('loading query kmers...', end=' ')
            bf = khmer.Nodetable(ksize, 1, 1)

            query_kmers = set()
            for record in screed.open(query):
                query_kmers.update(bf.get_kmer_hashes(record.sequence))

            print('got {}'.format(len(query_kmers)))

            # construct dict cdbg_id -> # of query k-mers
            cdbg_match_counts = kmer_idx.get_match_counts(query_kmers)
            for k, v in cdbg_match_counts.items():
                assert v <= kmer_idx.get_cdbg_size(k), k

            total_match_kmers = sum(cdbg_match_counts.values())
            f_found = total_match_kmers / len(query_kmers)
            print('=> containment: {:.1f}%'.format(f_found * 100))
            print('done loading & counting query k-mers in cDBG.'
                  ' ({:.1f}s)'.format(time.time() - start_time))

            total_kmers_in_cdbg_matches = 0
            for cdbg_id in set(cdbg_match_counts.keys()):
                total_kmers_in_cdbg_matches += kmer_idx.get_cdbg_size(cdbg_id)

            cdbg_sim = total_match_kmers / total_kmers_in_cdbg_matches
            print('cdbg match node similarity: {:.1f}%'.format(cdbg_sim * 100))
            cdbg_min_overhead = (total_kmers_in_cdbg_matches -
                                 total_match_kmers) /\
                total_kmers_in_cdbg_matches
            print('min cdbg overhead: {}'.format(cdbg_min_overhead))

            # calculate the cDBG matching k-mers sizes for each catlas node.
            catlas_match_counts =\
                kmer_idx.build_catlas_match_counts(cdbg_match_counts, dag,
                                                   dag_levels, layer1_to_cdbg)

            # check a few things - we've propogated properly:
            assert sum(cdbg_match_counts.values()) == \
                catlas_match_counts[top_node_id]
            # ...and all nodes have no more matches than total k-mers.
            for k, v in catlas_match_counts.items():
                assert v <= node_sizes[k], k

            # calculate the minimum overhead of the search, based on level 1
            # nodes.
            catlas_min_overhead = 0
            if catlas_match_counts[top_node_id]:
                all_query_kmers = catlas_match_counts[top_node_id]
                total_kmers_in_query_nodes = 0
                for node_id, level in dag_levels.items():
                    if level == 1 and catlas_match_counts.get(node_id):
                        total_kmers_in_query_nodes += node_sizes[node_id]

                catlas_min_overhead = (total_kmers_in_query_nodes -
                                       all_query_kmers) /\
                    total_kmers_in_query_nodes
                print('minimum catlas overhead: {}'.format(
                    catlas_min_overhead))

            # gather results of all queries
            fuzzy = args.max_overhead != 1.0
            if fuzzy:
                max_oh = args.max_overhead
                min_con = args.min_containment
                total_frontier = collect_frontier(dag, top_node_id,
                                                  node_sizes,
                                                  catlas_match_counts,
                                                  max_overhead=max_oh,
                                                  min_containment=min_con)
            else:
                total_frontier = collect_frontier_exact(dag, top_node_id,
                                                        node_sizes,
                                                        catlas_match_counts,
                                                        overhead=args.overhead,
                                                        verbose=args.verbose)

            # calculate level 1 nodes for this frontier in the catlas
            total_shadow = find_shadow(total_frontier, dag)

            # calculate associated cDBG nodes
            cdbg_shadow = set()
            for x in total_shadow:
                cdbg_shadow.update(layer1_to_cdbg.get(x))

            # done with main loop! now extract contigs using cDBG shadow
            # node list.
            print('done searching! {} frontier, {} catlas shadow nodes, {}'
                  ' cdbg nodes.'.format(len(total_frontier), len(total_shadow),
                                        len(cdbg_shadow)))

            # track extracted info
            total_bp = 0
            total_seqs = 0

            # build check MinHash w/seed=42
            query_sig = build_query_mh_for_seed(42, ksize, scaled, query)

            # track minhash of retrieved contigs using original query minhash:
            contigs_minhash = query_sig.minhash.copy_and_clear()

            retrieve_start = time.time()

            # walk through the contigs, retrieving.
            print('extracting contigs...')
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
            print('...fetched {} contigs, {} bp matching combined frontiers. '
                  ' ({:.1f}s)'.format(total_seqs, total_bp, time.time() -
                                      retrieve_start))

            # calculate summary values of extracted contigs
            containment = query_sig.minhash.contained_by(contigs_minhash)
            similarity = query_sig.minhash.similarity(contigs_minhash)
            print('query inclusion by retrieved contigs:'
                  ' {:.3f}%'.format(containment*100))
            print('query similarity to retrieved contigs:'
                  ' {:.3f}%'.format(similarity*100))

            num_seeds = 0

            # recover from above.
            best_containment = f_found

            # output to results.csv!
            csv_writer.writerow([query, containment, similarity, total_bp,
                                 total_seqs, ksize, len(query_kmers),
                                 best_containment, cdbg_min_overhead,
                                 catlas_min_overhead])
            csvoutfp.flush()

            # write out signature from retrieved contigs.
            sig_filename = os.path.basename(query) + '.contigs.sig'
            with open(os.path.join(args.output, sig_filename), 'wt') as fp:
                ss = sourmash_lib.SourmashSignature(contigs_minhash,
                                                    name=sig_filename)
                sourmash_lib.save_signatures([ss], fp)

            # write out cDBG IDs
            cdbg_listname = os.path.basename(query) + '.cdbg_ids.txt.gz'
            with gzip.open(os.path.join(args.output, cdbg_listname),
                           'wt') as fp:
                fp.write("\n".join([str(x) for x in cdbg_shadow]))

            # write out frontier nodes by seed
            frontier_listname = os.path.basename(query) + '.frontier.txt.gz'
            with gzip.open(os.path.join(args.output, frontier_listname),
                           'wt') as fp:
                for node, seedlist in total_frontier.items():
                    fp.write('{},{}\n'.format(node,
                                              " ".join([str(x) for x in
                                                        seedlist])))

            # write response curve
            response_curve_filename = os.path.basename(query) + '.response.txt'
            response_curve_filename = os.path.join(args.output,
                                                   response_curve_filename)
            search_utils.output_response_curve(response_curve_filename,
                                               cdbg_match_counts,
                                               kmer_idx,
                                               layer1_to_cdbg)
            print('total time: {:.1f}s'.format(time.time() - start_time))
        except KeyboardInterrupt:
            raise
        except:
            traceback.print_exc()

    # end main loop!

    sys.exit(0)


if __name__ == '__main__':
    main(sys.argv[1:])
