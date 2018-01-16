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

import screed

import sourmash_lib
from sourmash_lib import MinHash
from sourmash_lib.sourmash_args import load_query_signature
from sourmash_lib._minhash import hash_murmur

from .search_utils import (get_minhashdb_name, get_reads_by_cdbg,
                           build_queries_for_seeds, parse_seeds_arg)
from .search_utils import MinhashSqlDB
from spacegraphcats.logging import log
from search.frontier_search import (frontier_search_exact, find_shadow, NoContainment)
from . import search_utils
from .search_utils import (load_dag, load_layer1_to_cdbg, load_minhash, SqlKmerIndex)


def collect_frontier(dag, top_node_id,
                     node_kmer_sizes, node_query_kmers,
                     overhead=0.0, verbose=False):
    start = time.time()

    # gather results into total_frontier
    total_frontier = defaultdict(set)

    # do queries!
    try:
        frontier, num_leaves, num_empty, frontier_mh = \
          frontier_search_exact(top_node_id, dag, node_kmer_sizes, node_query_kmers, overhead)
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
    p.add_argument('-v', '--verbose', action='store_true')

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

    # load mapping between dom nodes and cDBG/graph nodes:
    layer1_to_cdbg = load_layer1_to_cdbg(catlas_to_cdbg, domfile)
    print('loaded {} layer 1 catlas nodes'.format(len(layer1_to_cdbg)))

    # find the contigs filename
    contigs = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

    # ...and index.
    kmer_sqlindex_file = contigs + '.sqlindex'
    kmer_idx = SqlKmerIndex(kmer_sqlindex_file)

    # load kmers per cdbg node
    print('loading sizes...')
    cdbg_kmer_sizes = {}
    for cdbg_id, size in kmer_idx.get_node_sizes():
        cdbg_kmer_sizes[cdbg_id] = size
    print('...done')

    # calculate the cDBG shadow sizes for each catlas node.
    x = []
    for (node_id, level) in dag_levels.items():
        x.append((level, node_id))
    x.sort()

    node_kmer_sizes = {}
    for level, node_id in x:
        if level == 1:
            total_kmers = 0
            for cdbg_node in layer1_to_cdbg.get(node_id):
                total_kmers += cdbg_kmer_sizes[node_id]
            
            node_kmer_sizes[node_id] = total_kmers
        else:
            sub_size = 0
            for child_id in dag[node_id]:
                sub_size += node_kmer_sizes[child_id]
            node_kmer_sizes[node_id] = sub_size

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
    csv_writer.writerow(['query', 'containment', 'similarity', 'bp', 'reads', 'n_seeds', 'ksize', 'scaled', 'best_containment'])

    # iterate over each query, do the thing.
    for query in args.query:
        # ignore all the problems!
        try:
            print('------')
            print('QUERY FILE:', query)
            start_time = time.time()

            # build a bf for the query
            print('loading query kmers...')
            bf = khmer.Nodetable(ksize, 1, 1)

            cdbg_count = defaultdict(int)

            x = set()
            n = 0
            for record in screed.open(query):
                for hashval in bf.get_kmer_hashes(record.sequence):
                    n += 1
                    if n % 250000 == 0:
                        print('...', n)
                    cdbg_id = kmer_idx.retrieve_sample_by_kmer(hashval)
                    cdbg_count[cdbg_id] += 1

            print('...done.')
            del cdbg_count[None]
            print('XXX', sum(cdbg_count.values()), len(cdbg_count))


            # calculate the cDBG matching k-mers sizes for each catlas node.
            x = []
            for (node_id, level) in dag_levels.items():
                x.append((level, node_id))
            x.sort()

            node_query_kmers = {}
            for level, node_id in x:
                if level == 1:
                    query_kmers = 0
                    for cdbg_node in layer1_to_cdbg.get(node_id):
                        query_kmers += cdbg_count.get(cdbg_node, 0)

                    if query_kmers:
                        node_query_kmers[node_id] = query_kmers
                else:
                    sub_size = 0
                    for child_id in dag[node_id]:
                        sub_size += node_query_kmers.get(child_id, 0)

                    if sub_size:
                        node_query_kmers[node_id] = sub_size

            print('ZZZ', node_query_kmers[node_id])

            assert sum(cdbg_count.values()) == node_query_kmers[node_id]
            print('XXX', sum(cdbg_count.values()), len(cdbg_count))
            print('...', n)

            # gather results of all queries across all seeds
            total_frontier = collect_frontier(dag, top_node_id,
                                              node_kmer_sizes,
                                              node_query_kmers,
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
            print('done searching! {} frontier, {} catlas shadow nodes, {} cdbg nodes.'.format(len(total_frontier), len(total_shadow), len(cdbg_shadow)))

            # track extracted info
            total_bp = 0
            total_seqs = 0

            # build check MinHash w/seed=42
            query_sig = build_queries_for_seeds([42], ksize, scaled, query)[0]

            # track minhash of retrieved contigs using original query minhash:
            contigs_minhash = query_sig.minhash.copy_and_clear()

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
            print('...fetched {} contigs, {} bp matching combined frontiers.'.format(total_seqs, total_bp))

            # calculate summary values of extracted contigs
            containment = query_sig.minhash.contained_by(contigs_minhash)
            similarity = query_sig.minhash.similarity(contigs_minhash)
            print('query inclusion by retrieved contigs:', containment)
            print('query similarity to retrieved contigs:', similarity)

            num_seeds = 0

            # in var query, we no longer have the ability to calculate this!
            best_containment = 0

            # output to results.csv!
            csv_writer.writerow([query, containment, similarity, total_bp, total_seqs, num_seeds, ksize, scaled, best_containment])
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

            print('total time: {:.1f}s'.format(time.time() - start_time))
        except KeyboardInterrupt:
            raise
        except:
            traceback.print_exc()

    ## end main loop!

    sys.exit(0)


if __name__ == '__main__':
    main()
