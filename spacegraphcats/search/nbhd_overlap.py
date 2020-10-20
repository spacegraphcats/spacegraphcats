#! /usr/bin/env python
"""

"""
import argparse
import csv
import gzip
import os
import sys
import time

import khmer
import screed
import sourmash_lib
from sourmash_lib import MinHash

from .frontier_search import (collect_frontier, collect_frontier_exact)
from . import search_utils
# from .search_utils import load_dag, load_kmer_index, load_layer1_to_cdbg
from .index import MPHF_KmerIndex
from .catlas import CAtlas


class QueryOutput:
    def __init__(self, query, catlas, kmer_idx, frontier, scaled):
        self.query = query
        self.catlas = catlas
        self.kmer_idx = kmer_idx
        self.frontier = frontier
        self.shadow = self.catlas.shadow(self.frontier)
        self.total_bp = 0
        self.total_seq = 0
        self.query_sig = self.query.build_query_mh_for_seed(42, scaled)
        self.contigs_minhash = self.query_sig.minhash.copy_and_clear()

    def __add_sequence(self, sequence):
        self.contigs_minhash.add_sequence(str(sequence), True)
        self.total_bp += len(sequence)
        self.total_seq += 1

    def containment(self):
        return self.query_sig.minhash.contained_by(self.contigs_minhash)

    def similarity(self):
        return self.query_sig.minhash.similarity(self.contigs_minhash)

    def retrieve_contigs(self, contigs):
        "extract contigs using cDBG shadow."

        # node list.
        # track extracted info
        retrieve_start = time.time()
        # walk through the contigs, retrieving.
        print('extracting contigs...')
        for n, record in enumerate(
                search_utils.get_contigs_by_cdbg(contigs,
                                                 self.shadow)):
            if n and n % 10000 == 0:
                offset_f = self.total_seq / len(self.shadow)
                print('...at n {} ({:.1f}% of shadow)'.format(self.total_seq,
                      offset_f * 100), end='\r')

            # track retrieved sequences in a minhash
            self.__add_sequence(record.sequence)

        # done - got all contigs!
        print('...fetched {} contigs, {} bp matching combined frontiers. '
              ' ({:.1f}s)'.format(self.total_seq, self.total_bp,
                                  time.time() - retrieve_start))

        # calculate summary values of extracted contigs
        print('query inclusion by retrieved contigs:'
              ' {:.3f}%'.format(self.containment()*100))
        print('query similarity to retrieved contigs:'
              ' {:.3f}%'.format(self.similarity()*100))

    def write(self, csv_writer, csvoutfp, outdir):
        containment = self.containment()
        similarity = self.similarity()
        q_name = self.query.filename
        bp = self.total_bp
        seqs = self.total_seq
        k = self.query.ksize
        num_q_kmers = len(self.query.kmers)
        (best_con,
         cdbg_min_oh,
         catlas_min_oh) = self.query.con_sim_upper_bounds(self.catlas,
                                                          self.kmer_idx)
        # output to results.csv!
        csv_writer.writerow([q_name, containment, similarity, bp,
                             seqs, k, num_q_kmers,
                             best_con, cdbg_min_oh,
                             catlas_min_oh])
        csvoutfp.flush()

        # write out signature from retrieved contigs.
        sig_filename = os.path.basename(q_name) + '.contigs.sig'
        with open(os.path.join(outdir, sig_filename), 'wt') as fp:
            ss = sourmash_lib.SourmashSignature(self.contigs_minhash,
                                                name='nbhd:'+self.query.name,
                                                filename=sig_filename)
            sourmash_lib.save_signatures([ss], fp)

        # write out cDBG IDs
        cdbg_listname = os.path.basename(q_name) + '.cdbg_ids.txt.gz'
        with gzip.open(os.path.join(outdir, cdbg_listname), 'wt') as fp:
            fp.write("\n".join([str(x) for x in sorted(self.shadow)]))

        # write out frontier nodes by seed
        frontier_listname = os.path.basename(q_name) + '.frontier.txt.gz'
        with gzip.open(os.path.join(outdir, frontier_listname), 'wt') as fp:
            for node, seedlist in sorted(self.frontier.items()):
                fp.write('{},{}\n'.format(node,
                                          " ".join([str(x) for x in
                                                    sorted(seedlist)])))

        # write response curve
        response_curve_filename = os.path.basename(q_name) + '.response.txt'
        response_curve_filename = os.path.join(outdir,
                                               response_curve_filename)
        cdbg_match_counts = self.query.cdbg_match_counts[self.catlas.name]
        search_utils.output_response_curve(response_curve_filename,
                                           cdbg_match_counts,
                                           self.kmer_idx,
                                           self.catlas.layer1_to_cdbg)


class Query:
    def __init__(self, query_file, ksize):
        self.filename = query_file
        self.ksize = ksize
        self.kmers = set()
        self.name = None
        print('----')
        print('QUERY FILE:', self.filename)

        # build hashes for all the query k-mers
        print('loading query kmers...', end=' ')
        bf = khmer.Nodetable(ksize, 1, 1)

        for record in screed.open(self.filename):
            if self.name is None:
                self.name = record.name
            self.kmers.update(bf.get_kmer_hashes(record.sequence))

        print('got {}'.format(len(self.kmers)))

        self.cdbg_match_counts = {}
        self.catlas_match_counts = {}

    def build_query_mh_for_seed(self, seed, scaled):
        mh = MinHash(0, self.ksize, scaled=scaled, seed=seed)
        for record in screed.open(self.filename):
            mh.add_sequence(record.sequence, True)

        return sourmash_lib.SourmashSignature(mh, name=self.name,
                                              filename=self.filename)

    def execute(self, catlas, kmer_idx, scaled, overhead, verbose,
                max_overhead, min_containment):
        cat_id = catlas.name
        # construct dict cdbg_id -> # of query k-mers
        self.cdbg_match_counts[cat_id] = kmer_idx.get_match_counts(self.kmers)
        for k, v in self.cdbg_match_counts[cat_id].items():
            assert v <= kmer_idx.get_cdbg_size(k), k
        # calculate the cDBG matching k-mers sizes for each catlas node.
        self.catlas_match_counts[cat_id] =\
            kmer_idx.build_catlas_match_counts(self.cdbg_match_counts[cat_id],
                                               catlas)

        # check a few things - we've propagated properly:
        assert sum(self.cdbg_match_counts[cat_id].values()) ==\
            self.catlas_match_counts[cat_id][catlas.root]
        # ...and all nodes have no more matches than total k-mers.
        for v, match_amount in self.catlas_match_counts[cat_id].items():
            assert match_amount <= catlas.index_sizes[v], v

        self.con_sim_upper_bounds(catlas, kmer_idx)

        # gather results of all queries
        fuzzy = max_overhead != 1.0
        if fuzzy:
            frontier = collect_frontier(catlas,
                                        self.catlas_match_counts[cat_id],
                                        max_overhead=max_overhead,
                                        min_containment=min_containment)
        else:
            frontier = collect_frontier_exact(catlas,
                                              self.catlas_match_counts[cat_id],
                                              overhead=overhead,
                                              verbose=verbose)

        # calculate level 1 nodes for this frontier in the catlas
        leaves = catlas.leaves(frontier)
        # calculate associated cDBG nodes
        cdbg_shadow = catlas.shadow(leaves)

        print('done searching! {} frontier, {} catlas shadow nodes, {}'
              ' cdbg nodes.'.format(len(frontier), len(leaves),
                                    len(cdbg_shadow)))
        return QueryOutput(self, catlas, kmer_idx, frontier, scaled)

    def con_sim_upper_bounds(self, catlas, kmer_idx):
        """
        Compute "best possible" bounds on the containment and similarity of the
        output of this query."""
        root = catlas.root

        cdbg_match_counts = self.cdbg_match_counts[catlas.name]
        catlas_match_counts = self.catlas_match_counts[catlas.name]
        total_match_kmers = sum(cdbg_match_counts.values())
        best_containment = total_match_kmers / len(self.kmers)
        print('=> containment: {:.1f}%'.format(best_containment * 100))

        total_kmers_in_cdbg_matches = 0
        for cdbg_id in set(cdbg_match_counts.keys()):
            total_kmers_in_cdbg_matches += kmer_idx.get_cdbg_size(cdbg_id)

        cdbg_sim = total_match_kmers / total_kmers_in_cdbg_matches
        print('cdbg match node similarity: {:.1f}%'.format(cdbg_sim * 100))
        cdbg_min_overhead = (total_kmers_in_cdbg_matches-total_match_kmers) /\
            total_kmers_in_cdbg_matches
        print('min cdbg overhead: {}'.format(cdbg_min_overhead))

        # calculate the minimum overhead of the search, based on level 1
        # nodes.
        catlas_min_overhead = 0
        if catlas_match_counts[root]:
            all_query_kmers = catlas_match_counts[root]
            total_kmers_in_query_nodes = 0
            for node_id, level in catlas.levels.items():
                if level == 1 and catlas_match_counts.get(node_id):
                    total_kmers_in_query_nodes += catlas.index_sizes[node_id]

            catlas_min_overhead = (total_kmers_in_query_nodes -
                                   all_query_kmers) /\
                total_kmers_in_query_nodes
            print('minimum catlas overhead: {}'.format(catlas_min_overhead))

        return best_containment, cdbg_min_overhead, catlas_min_overhead


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('query1')
    p.add_argument('query2')
    p.add_argument('-k', '--ksize', default=31, type=int,
                   help='k-mer size (default: 31)')
    p.add_argument('--scaled', default=1000, type=float,
                   help="scaled value for contigs minhash output")
    p.add_argument('-v', '--verbose', action='store_true')

    args = p.parse_args(argv)

    # figure out catlas and domfile information.
    catlas_file = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # load catlas DAG
    catlas = CAtlas(catlas_file, domfile=domfile)
    print('loaded {} nodes from catlas {}'.format(len(catlas), catlas_file))
    print('loaded {} layer 1 catlas nodes'.format(len(catlas.layer1_to_cdbg)))

    # find the contigs filename
    contigs = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

    # ...and kmer index.
    ki_start = time.time()
    kmer_idx = MPHF_KmerIndex.from_catlas_directory(args.catlas_prefix)
    print('loaded {} k-mers in index ({:.1f}s)'.format(
                                len(kmer_idx.mphf_to_kmer),
                                time.time() - ki_start))

    # calculate the k-mer sizes for each catlas node.
    catlas.decorate_with_index_sizes(kmer_idx)

    # get a single ksize & scaled
    ksize = int(args.ksize)
    scaled = int(args.scaled)

    # iterate over each query, do the thing.
    query1 = Query(args.query1, ksize)
    query2 = Query(args.query2, ksize)

    q1_output = query1.execute(catlas, kmer_idx, scaled, 0.0, args.verbose, 1.0, 1.0)
    q2_output = query2.execute(catlas, kmer_idx, scaled, 0.0, args.verbose, 1.0, 1.0)

    q1_output.retrieve_contigs(contigs)
    q2_output.retrieve_contigs(contigs)

    print('similarity', q1_output.contigs_minhash.similarity(q2_output.contigs_minhash))
    print('q1g cnt by q2', q1_output.query_sig.minhash.contained_by(q2_output.contigs_minhash))
    print('q2g cnt by q1', q2_output.query_sig.minhash.contained_by(q1_output.contigs_minhash))

    print('q1 cnt by q2', q1_output.query_sig.minhash.contained_by(q2_output.query_sig.minhash))

    sys.exit(0)


if __name__ == '__main__':
    main(sys.argv[1:])
