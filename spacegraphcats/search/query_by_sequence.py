#! /usr/bin/env python
"""
Query a catlas with a sequence (read, contig, or genome), and retrieve
cDBG node IDs and MinHash signatures for the matching unitigs in the graph.

XX optionally retrieve/save contig sequences? It wouldn't add much time...
Just disk space.
"""
import argparse
import csv
import gzip
import os
import sys
import time
import sqlite3

import screed
import sourmash
from sourmash import MinHash

from ..utils.logging import notify, error
from . import search_utils
from . import MPHF_KmerIndex, hash_sequence
from .catlas import CAtlas


class QueryOutput:
    def __init__(self, query, catlas, kmer_idx, doms, catlas_name):
        self.query = query
        self.catlas = catlas
        self.kmer_idx = kmer_idx
        self.doms = doms
        # here is where we go from level 1/dominators to cdbg nodes.
        self.shadow = self.catlas.shadow(doms)
        self.total_bp = 0
        self.total_seq = 0
        self.query_sig = self.query.sig
        self.contigs_minhash = self.query.mh.copy_and_clear()
        self.catlas_name = catlas_name

    def __add_sequence(self, sequence):
        self.contigs_minhash.add_sequence(str(sequence), True)
        self.total_bp += len(sequence)
        self.total_seq += 1

    def containment(self):
        return self.query_sig.minhash.contained_by(self.contigs_minhash)

    def similarity(self):
        return self.query_sig.minhash.similarity(self.contigs_minhash)

    def retrieve_contigs(self, contigs_db):
        "extract contigs using cDBG shadow."

        # node list.
        # track extracted info
        retrieve_start = time.time()
        # walk through the contigs, retrieving.
        notify("extracting contigs...")

        contigs_iter = search_utils.get_contigs_by_cdbg_sqlite(contigs_db, self.shadow)
        for n, record in enumerate(contigs_iter):
            if n and n % 10000 == 0:
                offset_f = self.total_seq / len(self.shadow)
                notify(
                    "...at n {} ({:.1f}% of shadow)",
                    self.total_seq,
                    offset_f * 100,
                    end="\r",
                )

            # track retrieved sequences in a minhash
            self.__add_sequence(record.sequence)

        # done - got all contigs!
        notify(
            "...fetched {} contigs, {} bp matching combined frontiers. " " ({:.1f}s)",
            self.total_seq,
            self.total_bp,
            time.time() - retrieve_start,
        )

        # calculate summary values of extracted contigs
        notify(
            "query inclusion by retrieved contigs:" " {:.1f}%", self.containment() * 100
        )
        notify(
            "query similarity to retrieved contigs:" " {:.1f}%", self.similarity() * 100
        )

    def write(self, csv_writer, csvoutfp, outdir, catlas_name):
        containment = self.containment()
        similarity = self.similarity()
        q_name = self.query.filename
        bp = self.total_bp
        seqs = self.total_seq
        k = self.query.ksize
        num_q_kmers = len(self.query.kmers)
        (best_con, cdbg_min_oh, catlas_min_oh) = self.query.con_sim_upper_bounds(
            self.catlas, self.kmer_idx
        )
        # output to results.csv!
        csv_writer.writerow(
            [
                q_name,
                containment,
                similarity,
                bp,
                seqs,
                k,
                num_q_kmers,
                best_con,
                cdbg_min_oh,
                catlas_min_oh,
                catlas_name,
            ]
        )
        csvoutfp.flush()

        # write out signature from retrieved contigs.
        sig_filename = os.path.basename(q_name) + ".contigs.sig"
        sig_name = "nbhd:{} from {}".format(self.query.name, self.catlas_name)
        with open(os.path.join(outdir, sig_filename), "wt") as fp:
            ss = sourmash.SourmashSignature(
                self.contigs_minhash, name=sig_name, filename=sig_filename
            )
            sourmash.save_signatures([ss], fp)

        # write out cDBG IDs
        cdbg_listname = os.path.basename(q_name) + ".cdbg_ids.txt.gz"
        with gzip.open(os.path.join(outdir, cdbg_listname), "wt") as fp:
            id_list = "\n".join([str(x) for x in sorted(self.shadow)])
            print(id_list, file=fp)

        # write out catlas nodes
        frontier_listname = os.path.basename(q_name) + ".frontier.txt.gz"
        with gzip.open(os.path.join(outdir, frontier_listname), "wt") as fp:
            for node in sorted(self.doms):
                fp.write("{}\n".format(node))

        # write response curve
        response_curve_filename = os.path.basename(q_name) + ".response.txt"
        response_curve_filename = os.path.join(outdir, response_curve_filename)
        cdbg_match_counts = self.query.cdbg_match_counts
        search_utils.output_response_curve(
            response_curve_filename,
            cdbg_match_counts,
            self.kmer_idx,
            self.catlas.layer1_to_cdbg,
        )


class Query:
    def __init__(self, query_file, ksize, scaled, catlas_name, debug=True):
        self.filename = query_file
        self.ksize = ksize
        self.kmers = set()
        self.name = None
        mh = MinHash(0, ksize, scaled=scaled)
        self.mh = mh
        self.catlas_name = catlas_name
        self.debug = debug

        notify("----")
        notify("QUERY FILE: {}", self.filename)

        # build hashes for all the query k-mers & create signature
        notify("loading query kmers...", end=" ")

        for record in screed.open(self.filename):
            if self.name is None:
                self.name = record.name
            if len(record.sequence) >= int(ksize):
                self.kmers.update(hash_sequence(record.sequence, self.ksize))
            mh.add_sequence(record.sequence, True)

        self.sig = sourmash.SourmashSignature(
            mh, name=self.name, filename=self.filename
        )

        notify("got {} k-mers from query", len(self.kmers))

        self.cdbg_match_counts = None
        self.catlas_match_counts = None

    def execute(self, catlas, kmer_idx):
        # construct dict cdbg_id -> # of query k-mers
        self.cdbg_match_counts = kmer_idx.count_cdbg_matches(self.kmers)
        # (confirm that fewer k-mers are returned than the node size)
        for k, v in self.cdbg_match_counts.items():
            assert v <= kmer_idx.get_cdbg_size(k), k

        # propogate the match counts to the dominators & the catlas
        self.catlas_match_counts = kmer_idx.count_catlas_matches(
            self.cdbg_match_counts, catlas
        )

        if self.debug and self.catlas_match_counts:
            # check a few things - we've propagated properly:
            assert (
                sum(self.cdbg_match_counts.values())
                == self.catlas_match_counts[catlas.root]
            )

            # ...and all nodes have no more matches than total k-mers.
            for v, match_amount in self.catlas_match_counts.items():
                assert match_amount <= catlas.index_sizes[v], v

            # calculate best possible.
            self.con_sim_upper_bounds(catlas, kmer_idx)

        # gather domset nodes matching query k-mers
        cdbg_nodes = set()
        doms = set()
        for cdbg_node in self.cdbg_match_counts:
            # nodes that matched in cDBG
            cdbg_nodes.add(cdbg_node)
            # dominator for this cDBG node
            doms.add(catlas.cdbg_to_layer1[cdbg_node])

        notify(
            "done searching! {} dominator nodes, {} cdbg nodes.",
            len(doms),
            len(cdbg_nodes),
        )
        return QueryOutput(self, catlas, kmer_idx, doms, self.catlas_name)

    def con_sim_upper_bounds(self, catlas, kmer_idx):
        """
        Compute "best possible" bounds on the containment and similarity of the
        output of this query."""
        root = catlas.root

        cdbg_match_counts = self.cdbg_match_counts
        catlas_match_counts = self.catlas_match_counts
        total_match_kmers = sum(cdbg_match_counts.values())
        best_containment = total_match_kmers / len(self.kmers)
        notify("=> containment: {:.1f}%", best_containment * 100)

        total_kmers_in_cdbg_matches = 0
        for cdbg_id in set(cdbg_match_counts.keys()):
            total_kmers_in_cdbg_matches += kmer_idx.get_cdbg_size(cdbg_id)

        cdbg_sim = 0
        cdbg_min_overhead = 0
        if total_kmers_in_cdbg_matches:
            cdbg_sim = total_match_kmers / total_kmers_in_cdbg_matches
            notify("cdbg match node similarity: {:.1f}%", cdbg_sim * 100)
            cdbg_min_overhead = (
                total_kmers_in_cdbg_matches - total_match_kmers
            ) / total_kmers_in_cdbg_matches
            notify("min cdbg overhead: {:.1f}%", cdbg_min_overhead * 100)

        # calculate the minimum overhead of the search, based on level 1
        # nodes.
        catlas_min_overhead = 0
        if catlas_match_counts and catlas_match_counts[root]:
            all_query_kmers = catlas_match_counts[root]
            total_kmers_in_query_nodes = 0
            for node_id, level in catlas.levels.items():
                if level == 1 and catlas_match_counts.get(node_id):
                    total_kmers_in_query_nodes += catlas.index_sizes[node_id]

            catlas_min_overhead = (
                total_kmers_in_query_nodes - all_query_kmers
            ) / total_kmers_in_query_nodes
            notify("minimum catlas overhead: {:.1f}%", catlas_min_overhead * 100)

        return best_containment, cdbg_min_overhead, catlas_min_overhead


def main(argv):
    """\
    Query a catlas with a sequence (read, contig, or genome), and retrieve
    cDBG node IDs and MinHash signatures for the matching unitigs in the graph.
    """

    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("output")
    p.add_argument("--query", help="query sequences", nargs="+")
    p.add_argument("--contigs-db", required=True)
    p.add_argument(
        "-k", "--ksize", default=31, type=int, help="k-mer size (default: 31)"
    )
    p.add_argument(
        "--scaled",
        default=1000,
        type=float,
        help="scaled value for contigs minhash output",
    )
    p.add_argument("-v", "--verbose", action="store_true")
    #    p.add_argument('--cdbg-only', action='store_true',
    #                   help="(for paper eval) do not expand query using domset)")

    args = p.parse_args(argv)
    outdir = args.output

    if not args.query:
        print("must specify at least one query file using --query.")
        sys.exit(-1)

    # make sure all of the query sequences exist.
    for filename in args.query:
        if not os.path.exists(filename):
            error("query seq file {} does not exist.", filename)
            sys.exit(-1)

    # create output directory if it doesn't exist.
    try:
        os.mkdir(outdir)
    except OSError:
        pass
    if not os.path.isdir(outdir):
        error("output {} is not a directory", outdir)
        sys.exit(-1)

    # load catlas DAG
    catlas = CAtlas(args.catlas_prefix)
    notify("loaded {} nodes from catlas {}", len(catlas), args.catlas_prefix)
    notify("loaded {} layer 1 catlas nodes", len(catlas.layer1_to_cdbg))
    catlas_name = os.path.basename(args.catlas_prefix.rstrip(("/")))

    # ...and kmer index.
    ki_start = time.time()
    kmer_idx = MPHF_KmerIndex.from_catlas_directory(args.catlas_prefix)
    notify("loaded {} k-mers in index ({:.1f}s)", len(kmer_idx), time.time() - ki_start)

    # ...and contigs db
    contigs_db = sqlite3.connect(args.contigs_db)

    # calculate the k-mer sizes for each catlas node.
    catlas.decorate_with_index_sizes(kmer_idx)

    # get a single ksize & scaled
    ksize = int(args.ksize)
    scaled = int(args.scaled)

    # record command line
    with open(os.path.join(outdir, "command.txt"), "wt") as fp:
        fp.write(str(sys.argv))
        fp.write("\n")

    # output results.csv in the output directory:
    csvoutfp = open(os.path.join(outdir, "results.csv"), "wt")
    csv_writer = csv.writer(csvoutfp)
    csv_writer.writerow(
        [
            "query",
            "containment",
            "similarity",
            "bp",
            "contigs",
            "ksize",
            "num_query_kmers",
            "best_containment",
            "cdbg_min_overhead",
            "catlas_min_overhead",
            "catlas_name",
        ]
    )

    # iterate over each query, do the thing.
    for query_file in args.query:
        start_time = time.time()
        query = Query(query_file, ksize, scaled, catlas_name)
        if not len(query.kmers):
            notify("query {} is empty; skipping.", query_file)
            continue

        q_output = query.execute(catlas, kmer_idx)
        q_output.retrieve_contigs(contigs_db)
        notify("total time: {:.1f}s", time.time() - start_time)

        q_output.write(csv_writer, csvoutfp, outdir, args.catlas_prefix)
    # end main loop!

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
