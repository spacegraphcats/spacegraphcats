#! /usr/bin/env python
"""
Query the cdbg with a sequence (read, contig, or genome), and retrieve
contigs and/or reads for the matching unitigs in the graph.

Need to have run `index_cdbg` and maybe `index_reads` rules.
"""
import argparse
import os
import sys
import time
import pickle
from collections import defaultdict
import sqlite3

import screed

from . import search_utils
from ..utils.logging import notify, error
from . import MPHF_KmerIndex, hash_sequence
from .catlas import CAtlas


def main(argv):
    """\
    Query the cdbg with a sequence (read, contig, or genome), and retrieve
    contigs and/or reads for the matching unitigs in the graph.
    """

    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument("base_prefix", help="base dir")
    p.add_argument("cdbg_prefix", help="cdbg dir")
    p.add_argument("catlas_prefix", help="catlas dir")
    p.add_argument("query", help="query sequences", nargs="+")
    p.add_argument(
        "-k", "--ksize", default=31, type=int, help="k-mer size (default: 31)"
    )
    p.add_argument('-R', '--output-reads', action='store_true',
                   help="Output raw reads instead of cDBG unitigs")
    p.add_argument('-r', '--radius', type=int, default=0,
                   help="expand found nodes by this radius")
    p.add_argument("-v", "--verbose", action="store_true")

    args = p.parse_args(argv)

    retrieve_contigs = not args.output_reads

    # make sure all of the query sequences exist.
    for filename in args.query:
        if not os.path.exists(filename):
            error("query seq file {} does not exist.", filename)
            sys.exit(-1)

    if retrieve_contigs:
        contigs_db_file = os.path.join(args.cdbg_prefix, 'bcalm.unitigs.db')
        assert os.path.exists(contigs_db_file), f'{contigs_db_file} must exist'
    else:                                 # retrieve reads
        reads_idx_file = os.path.join(args.catlas_prefix, 'reads.bgz.index')
        assert os.path.exists(reads_idx_file), f'{reads_idx_file} must exist'

        reads_bgz_file = os.path.join(args.base_prefix, f'{args.base_prefix}.reads.bgz')
        assert os.path.exists(reads_bgz_file), f'{reads_bgz_file} must exist'

    # ...and kmer index.
    ki_start = time.time()
    kmer_idx = MPHF_KmerIndex.from_catlas_directory(args.catlas_prefix)
    notify("loaded {} k-mers in index ({:.1f}s)", len(kmer_idx), time.time() - ki_start)

    # get ksize
    ksize = int(args.ksize)

    matching_cdbg = set()
    for filename in args.query:
        notify(f"Reading from query '{filename}'")
        for record in screed.open(filename):
            if len(record.sequence) < int(ksize):
                continue

            kmers = hash_sequence(record.sequence, ksize)
            cdbg_match_counts = kmer_idx.count_cdbg_matches(kmers)

            notify(
                f"got {len(cdbg_match_counts)} matching cdbg nodes for '{record.name[:15]}...' ({len(kmers)} kmers)"
            )

            matching_cdbg.update(cdbg_match_counts)

    notify(f"Found {len(matching_cdbg)} matching cDBG nodes total for all query sequences.")

    if args.radius:
        notify(f"inflating by radius {args.radius}...")

        gxtfile = os.path.join(args.catlas_prefix, "cdbg.gxt")
        graph_neighbors = defaultdict(set)
        with open(gxtfile, 'rt') as fp:
            fp.readline()                 # discard first
            for edgeline in fp:
                a, b = edgeline.strip().split()
                a = int(a)
                b = int(b)

                graph_neighbors[a].add(b)
                graph_neighbors[b].add(a)

        orig_matching_cdbg = matching_cdbg

        for i in range(args.radius):
            notify(f"round {i+1} inflation...")
            new_nodes = set(matching_cdbg)
            for node in matching_cdbg:
                new_nodes.update(graph_neighbors[node])

            notify(f"...inflated from {len(matching_cdbg)} to {len(new_nodes)} nodes.")
            matching_cdbg = new_nodes

    if retrieve_contigs:
        notify("Retrieving contigs...")
        contigs_db = sqlite3.connect(contigs_db_file)
        contigs_iter = search_utils.get_contigs_by_cdbg_sqlite(contigs_db,
                                                               matching_cdbg)

        n = -1
        for n, record in enumerate(contigs_iter):
            print(f">{record.name}\n{record.sequence}")

        notify(f"...output {n+1} contigs to stdout.")
    else:
        notify("Retrieving reads...")
        reads_index = sqlite3.connect(reads_idx_file)
        reads_iter = search_utils.get_reads_by_cdbg(reads_idx_file,
                                                    reads_bgz_file,
                                                    matching_cdbg)

        n = -1
        for n, (record, offset) in enumerate(reads_iter):
            print(f">{record.name}\n{record.sequence}")

        notify(f"...output {n+1} reads to stdout.")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
