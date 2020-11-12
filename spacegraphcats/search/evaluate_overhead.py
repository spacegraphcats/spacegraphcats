#! /usr/bin/env python
"""
Evaluate the minimum overhead of a query into a catlas.
"""
import argparse
import os
import sys
import gzip
import sqlite3

import khmer

import screed

from .catlas import CAtlas
from spacegraphcats.search import search_utils


def main(argv):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("catlas_prefix")
    p.add_argument("query")
    p.add_argument("cdbg_nodefile")
    p.add_argument("--contigs-db", required=True)

    p.add_argument("-o", "--output", type=argparse.FileType("wt"))
    p.add_argument(
        "-k", "--ksize", default=31, type=int, help="k-mer size (default: 31)"
    )
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    contigs_db = sqlite3.connect(args.contigs_db)

    assert args.output, "must specify -o"
    outfp = args.output
    outname = args.output.name

    print("loading bf...", end=" ")
    bf = khmer.Nodetable(args.ksize, 3e8, 2)
    bf.consume_seqfile(args.query)
    print("done.")

    print("loading catlas...", end=" ")
    catlas = CAtlas(args.catlas_prefix)
    layer1_to_cdbg = catlas.layer1_to_cdbg
    print("done.")

    print("loading nodefile {}".format(args.cdbg_nodefile))
    cdbg_nodes = set()
    with gzip.open(args.cdbg_nodefile, "r") as fp:
        for line in fp:
            cdbg_nodes.add(int(line.strip()))

    print("loading contigs")
    total_bp = 0
    total_seqs = 0

    n_homogeneous = 0
    n_missing = 0
    bp_missing = 0
    for n, record in enumerate(
        search_utils.get_contigs_by_cdbg_sqlite(contigs_db, cdbg_nodes)
    ):
        counts = bf.get_kmer_counts(record.sequence)
        if min(counts) == max(counts):
            n_homogeneous += 1

        if max(counts) == 0:
            n_missing += 1
            bp_missing += len(record.sequence)

        outfp.write("{}\n".format(len(record.sequence)))

        total_bp += len(record.sequence)
        total_seqs += 1

    print("")
    print("fetched {} contigs, {} bp matching node list.".format(total_seqs, total_bp))
    print("n_homogeneous: {}".format(n_homogeneous))
    print("pure overhead count: {} seqs / {} bp".format(n_missing, bp_missing))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
