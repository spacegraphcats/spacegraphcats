#! /usr/bin/env python
"""
Query the cdbg with a sequence (read, contig, or genome), and retrieve
contigs and/or reads for the matching unitigs in the graph.

'graphgrep_iter' does not need the k-mer index, unlike 'graphgrep'.

In order to retrieve reads, `index_reads` needs to have been run.
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

    # get ksize
    ksize = int(args.ksize)

    query_kmers = set()
    for filename in args.query:
        notify(f"Reading from query '{filename}'")
        for record in screed.open(filename):
            if len(record.sequence) < int(ksize):
                continue

            kmers = hash_sequence(record.sequence, ksize)

            notify(
                f"got {len(kmers)} kmers for '{record.name[:15]}...'"
            )
            query_kmers.update(kmers)

    notify(f"Loaded {len(query_kmers)} query kmers total.")

    # now do search to find primary nodes -

    matching_cdbg = set()

    unitigs_db = os.path.join(args.cdbg_prefix, 'bcalm.unitigs.db')
    db = sqlite3.connect(unitigs_db)
    for n, record in enumerate(search_utils.contigs_iter_sqlite(db)):
        if n and n % 1000 == 0:
            print(f'... {n} {len(matching_cdbg)} (unitigs)', end='\r', file=sys.stderr)

        record_kmers = set(hash_sequence(record.sequence, ksize))
        if record_kmers & query_kmers:
            matching_cdbg.add(record.name)

    print(f'loaded {n} query sequences (unitigs)', file=sys.stderr)
    print(f'total matches: {len(matching_cdbg)}', file=sys.stderr)

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
