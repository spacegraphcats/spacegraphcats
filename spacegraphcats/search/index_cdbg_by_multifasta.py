#! /usr/bin/env python
import argparse
import os
import sys
import time
import pickle
from collections import defaultdict

import screed

from ..utils.logging import notify, error
from . import MPHF_KmerIndex, hash_sequence
from .catlas import CAtlas


def main(argv):
    """\
    Query a catlas with a sequence (read, contig, or genome), and retrieve
    cDBG node IDs and MinHash signatures for the matching unitigs in the graph.
    """

    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("output")
    p.add_argument("--query", help="query sequences", nargs="+")
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

    args = p.parse_args(argv)
    outfile = args.output

    if not args.query:
        print("must specify at least one query file using --query.")
        sys.exit(-1)

    # make sure all of the query sequences exist.
    for filename in args.query:
        if not os.path.exists(filename):
            error("query seq file {} does not exist.", filename)
            sys.exit(-1)

    # load catlas DAG
    catlas = CAtlas(args.catlas_prefix)
    notify("loaded {} nodes from catlas {}", len(catlas), args.catlas_prefix)
    notify("loaded {} layer 1 catlas nodes", len(catlas.layer1_to_cdbg))

    # ...and kmer index.
    ki_start = time.time()
    kmer_idx = MPHF_KmerIndex.from_catlas_directory(args.catlas_prefix)
    notify("loaded {} k-mers in index ({:.1f}s)", len(kmer_idx), time.time() - ki_start)

    # calculate the k-mer sizes for each catlas node.
    catlas.decorate_with_index_sizes(kmer_idx)

    # get a single ksize & scaled
    ksize = int(args.ksize)
    scaled = int(args.scaled)

    records_to_cdbg = {}
    cdbg_to_records = defaultdict(set)
    for filename in args.query:
        print(f"Reading from '{filename}'")
        for record in screed.open(filename):
            if len(record.sequence) < int(ksize):
                continue

            kmers = hash_sequence(record.sequence, ksize)
            cdbg_match_counts = kmer_idx.count_cdbg_matches(kmers)

            print(
                f"got {len(cdbg_match_counts)} cdbg nodes for {record.name[:15]} ({len(kmers)} kmers)"
            )

            dominators = set()
            for cdbg_node in cdbg_match_counts:
                dominators.add(catlas.cdbg_to_layer1[cdbg_node])

            print(f"got {len(dominators)} dominators for {record.name[:15]}")

            shadow = catlas.shadow(dominators)
            print(f"got {len(shadow)} cdbg_nodes under {len(dominators)} dominators")

            records_to_cdbg[(filename, record.name)] = shadow
            for cdbg_node in shadow:
                cdbg_to_records[cdbg_node].add((filename, record.name))

    with open(outfile, "wb") as fp:
        print(f"saving pickled index to '{outfile}'")
        pickle.dump((args.catlas_prefix, records_to_cdbg, cdbg_to_records), fp)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
