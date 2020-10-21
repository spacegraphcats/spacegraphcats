#! /usr/bin/env python
"""
Benchmark the search algorithm, without I/O considerations.
"""
import sys
import os

# add spacegraphcats package to import path:
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import spacegraphcats
from spacegraphcats.search import CAtlas, MPHF_KmerIndex, Query
from collections import defaultdict
import argparse
import sys
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("project", help="Project directory", type=str)
    parser.add_argument("query", nargs="+")
    parser.add_argument("-o", "--output", type=str)
    parser.add_argument("-k", "--ksize", default=31, type=int)
    parser.add_argument("--validate", action="store_true")
    args = parser.parse_args()

    # figure out catlas and domfile information.
    catlas_file = os.path.join(args.project, "catlas.csv")
    domfile = os.path.join(args.project, "first_doms.txt")

    # load catlas DAG
    catlas = CAtlas(catlas_file, domfile=domfile)
    print("loaded {} nodes from catlas {}".format(len(catlas), catlas_file))
    print("loaded {} layer 1 catlas nodes".format(len(catlas.layer1_to_cdbg)))

    # ...and kmer index.
    ki_start = time.time()
    kmer_idx = MPHF_KmerIndex.from_catlas_directory(args.project)
    print(
        "loaded {} k-mers in index ({:.1f}s)".format(
            len(kmer_idx.mphf_to_kmer), time.time() - ki_start
        )
    )

    # calculate the k-mer sizes for each catlas node.
    catlas.decorate_with_index_sizes(kmer_idx)

    print("building index 2 (cdbg node ID to pieces) - untimed")
    cdbg_to_pieces = defaultdict(set)
    for node_id in catlas:
        level = catlas.levels[node_id]
        if level == 1:
            pieces = catlas.layer1_to_cdbg.get(node_id)
            for cdbg_node in pieces:
                cdbg_to_pieces[cdbg_node] = set(pieces)

    # get a single ksize
    ksize = int(args.ksize)

    # BENCHMARK:
    # iterate over each query, do the thing.
    total_kmers = 0
    total_matched_nodes = 0

    start = time.time()

    for query_file in args.query:
        query = Query(query_file, ksize)
        assert len(query.kmers)

        pieces = set()
        for kmer in query.kmers:
            cdbg_node = kmer_idx.get_cdbg_id(kmer)
            if cdbg_node is not None:
                pieces.update(cdbg_to_pieces[cdbg_node])

        print("found {} pieces for query.".format(len(pieces)))

        total_kmers += len(query.kmers)
        total_matched_nodes += len(pieces)

        if args.validate:
            print("validating...")
            q_output = query.execute(catlas, kmer_idx, 1000, 0.0, True, 1.0, 1.0)
            assert len(q_output.shadow) == total_matched_nodes, (
                q_output.shadow,
                total_matched_nodes,
            )
    # end main loop!

    end = time.time()

    # done. output.

    outfp = sys.stdout
    if args.output:
        outfp = open(args.output, "at")
    print(
        "{},{},{:.1f},{},search".format(
            total_kmers / len(args.query),
            total_matched_nodes / len(args.query),
            end - start,
            args.project,
        ),
        file=outfp,
    )


if __name__ == "__main__":
    sys.exit(main())
