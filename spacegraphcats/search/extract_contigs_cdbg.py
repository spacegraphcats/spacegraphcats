#! /usr/bin/env python
"""
Retrieve the unitigs containing a set of k-mers, using exact/MPHF code.

This _does not_ use the domset or neighborhood, merely recovers exactly the
unitigs that contain the query k-mers.
"""
import argparse
import os
import sys
import sqlite3

import screed

from spacegraphcats.cdbg import hash_sequence

from . import search_utils


def main(argv=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("query")
    p.add_argument("--contigs-db", required=True)
    p.add_argument(
        "-k", "--ksize", default=31, type=int, help="k-mer size (default: 31)"
    )
    p.add_argument("-o", "--output", type=argparse.FileType("wt"))
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    contigs_db = sqlite3.connect(args.contigs_db)

    # load k-mer MPHF index
    kmer_idx = search_utils.load_kmer_index(args.catlas_prefix)

    # build hashes for all the query k-mers
    print("loading query kmers...")

    x = set()
    n = 0

    query_kmers = set()
    for record in screed.open(args.query):
        query_kmers.update(hash_sequence(record.sequence, args.ksize))

    query_kmers = list(query_kmers)

    # find the list of cDBG nodes that contain at least one query k-mer
    cdbg_match_counts = kmer_idx.count_cdbg_matches(query_kmers)

    if not cdbg_match_counts:
        print("No k-mer matches!? (Check that you're using the right ksize.")
        return -1

    # calculate number of nodes found -
    cdbg_shadow = set(cdbg_match_counts.keys())

    # calculate the sum total k-mers across all of the matching nodes
    cdbg_node_sizes = {}
    for cdbg_id in cdbg_shadow:
        cdbg_node_sizes[cdbg_id] = kmer_idx.get_cdbg_size(cdbg_id)

    # output some stats
    total_found = sum(cdbg_match_counts.values())
    f_found = total_found / len(query_kmers)
    print("...done loading & counting query k-mers in cDBG.")
    print("containment: {:.1f}%".format(f_found * 100))

    total_kmers_in_cdbg_nodes = sum(cdbg_node_sizes.values())
    sim = total_found / total_kmers_in_cdbg_nodes
    print("similarity: {:.1f}%".format(sim * 100))

    if not args.output:
        sys.exit(0)

    # if output requested, extract unitigs.
    outfp = args.output
    outname = args.output.name

    total_bp = 0
    total_seqs = 0

    print("extracting contigs to {}.".format(outname))
    for n, record in enumerate(
        search_utils.get_contigs_by_cdbg_sqlite(contigs_db, cdbg_shadow)
    ):
        outfp.write(">{}\n{}\n".format(record.name, record.sequence))

        total_bp += len(record.sequence)
        total_seqs += 1

    print("")
    print("fetched {} contigs, {} bp matching node list.".format(total_seqs, total_bp))

    return 0


if __name__ == "__main__":
    sys.exit(main())
