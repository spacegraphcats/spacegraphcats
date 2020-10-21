#! /usr/bin/env python
"""
Benchmark the indexPieces algorithm, without I/O considerations.
"""
import sys
import os

# add spacegraphcats package to import path:
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import spacegraphcats
from spacegraphcats.search.catlas import CAtlas
from spacegraphcats.cdbg.index_contigs_by_kmer import build_mphf
from collections import defaultdict
import argparse
import sys
import time
import khmer
import screed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("project", help="Project directory", type=str)
    parser.add_argument("-o", "--output", type=str)
    parser.add_argument("-k", "--ksize", default=31, type=int)
    args = parser.parse_args()

    # figure out catlas and domfile information.
    catlas_file = os.path.join(args.project, "catlas.csv")
    domfile = os.path.join(args.project, "first_doms.txt")

    # grab contigs
    print("loading contigs...")
    contigs_filename = os.path.join(args.project, "contigs.fa.gz")
    contigs = {}
    for record in screed.open(contigs_filename):
        contigs[record.name] = record
    print("done. {} contigs.".format(len(contigs)))

    # load catlas DAG
    catlas = CAtlas(catlas_file, domfile=domfile)
    print("loaded {} nodes from catlas {}".format(len(catlas), catlas_file))
    print("loaded {} layer 1 catlas nodes".format(len(catlas.layer1_to_cdbg)))

    # create k-mer hashing machinery
    kh = khmer.Nodetable(args.ksize, 1, 1)

    # function to yield records as if reading from a file
    def create_records_iter():
        return contigs.values()

    # BENCHMARK:
    # build the two indices (kmer to cdbg node, cdbg node ID to pieces)
    start = time.time()
    print("building MPHF index.")
    x, mphf_to_kmer, mphf_to_cdbg, sizes = build_mphf(kh, create_records_iter)
    print("done! {:.1f}s.".format(time.time() - start))

    print("building index 2 (cdbg node ID to pieces)")
    cdbg_to_pieces = defaultdict(set)
    for node_id in catlas:
        level = catlas.levels[node_id]
        if level == 1:
            pieces = catlas.layer1_to_cdbg.get(node_id)
            for cdbg_node in pieces:
                cdbg_to_pieces[cdbg_node] = set(pieces)
    end = time.time()

    # done. output.

    outfp = sys.stdout
    if args.output:
        outfp = open(args.output, "at")
    print(
        "{},{},{:.1f},{},indexPieces".format(
            len(mphf_to_kmer), len(catlas), end - start, args.project
        ),
        file=outfp,
    )


if __name__ == "__main__":
    sys.exit(main())
