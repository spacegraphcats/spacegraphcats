#! /usr/bin/env python
"""
Use the cDBG and the info.csv file to estimate the abundance of a query
genome.

The query genome should be largely or completely enclosed in the ???
"""
import argparse
import os
import sys

import screed

from spacegraphcats.cdbg import hash_sequence
from . import search_utils


def main(argv=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("queries", nargs="+")
    p.add_argument(
        "-k", "--ksize", default=31, type=int, help="k-mer size (default: 31)"
    )
    p.add_argument("-o", "--output", type=argparse.FileType("wt"))
    args = p.parse_args(argv)

    assert args.output, "must supply -o"

    x = search_utils.load_cdbg_size_info(args.catlas_prefix)
    cdbg_kmer_sizes, cdbg_weighted_kmer_sizes = x

    # load k-mer MPHF index
    kmer_idx = search_utils.load_kmer_index(args.catlas_prefix)

    # build hashes for all the query k-mers
    print("loading query kmers...")

    print("queryfile,containment,mean_abundance", file=args.output)

    for query in args.queries:
        print("loading", query)
        query_kmers = set()
        for record in screed.open(query):
            hashes = hash_sequence(record.sequence, args.ksize)
            query_kmers.update(hashes)

        # find the list of cDBG nodes that contain at least one query k-mer
        cdbg_match_counts = kmer_idx.count_cdbg_matches(query_kmers)

        # calculate number of nodes found -
        cdbg_shadow = set(cdbg_match_counts.keys())

        # calculate the sum total k-mers across all of the matching nodes
        cdbg_node_sizes = {}
        cdbg_total_weighted = 0.0
        for cdbg_id in cdbg_shadow:
            cdbg_node_sizes[cdbg_id] = kmer_idx.get_cdbg_size(cdbg_id)
            cdbg_total_weighted += cdbg_weighted_kmer_sizes[cdbg_id]

        # output some stats
        total_found = sum(cdbg_match_counts.values())
        f_found = total_found / len(query_kmers)
        print("...done loading & counting query k-mers in cDBG.")
        print("containment: {:.1f}%".format(f_found * 100))

        weight = cdbg_total_weighted / total_found
        print("weight:", weight)

        if f_found < 0.5:
            print("skipping output for {}; low containment.".format(query))
            continue

        print("{},{},{}".format(query, f_found, weight), file=args.output)

    return 0


if __name__ == "__main__":
    sys.exit(main())
