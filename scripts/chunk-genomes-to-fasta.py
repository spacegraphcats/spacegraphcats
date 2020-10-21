#! /usr/bin/env python
import screed
import argparse
import numpy
from itertools import product
from sourmash._minhash import hash_murmur
from pickle import dump
import pickle

SIZE = 10000
OVERLAP = 100


def make_all(ksize):
    return ["".join(kmer) for kmer in product("ACGT", repeat=ksize)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes", nargs="+")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()

    assert args.output, "please specify -o"

    outfp = open(args.output, "wt")

    n = 0
    genome_n = 0
    group_info = {}
    group_ident = {}
    labels = {}
    node_id_to_group_idx = {}
    for genome in args.genomes:
        print(genome)
        genome_n += 1
        for record in screed.open(genome):
            for start in range(0, len(record.sequence), SIZE):
                outfp.write(
                    ">chunk{}\n{}\n.".format(
                        n, record.sequence[start : start + SIZE + OVERLAP]
                    )
                )
                n += 1


if __name__ == "__main__":
    main()
