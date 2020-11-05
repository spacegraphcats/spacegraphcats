#! /usr/bin/env python
import screed
import argparse
import sourmash
import numpy
from itertools import product
from sourmash._minhash import hash_murmur
from pickle import dump
import pickle

SIZE = 10000


def make_all(ksize):
    return ["".join(kmer) for kmer in product("ACGT", repeat=ksize)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genomes", nargs="+")
    parser.add_argument("-o", "--output")
    parser.add_argument(
        "-k", "--ksize", default=5, type=int, help="k-mer size for vectors"
    )
    args = parser.parse_args()

    assert args.output, "please specify -o"

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
                mh = sourmash.MinHash(
                    n=0, ksize=args.ksize, scaled=1, track_abundance=1
                )
                mh.add_sequence(record.sequence[start : start + SIZE], True)
                group_info[n] = mh

                mh = sourmash.MinHash(n=0, ksize=31, scaled=1000)
                mh.add_sequence(record.sequence[start : start + SIZE], True)
                group_ident[n] = mh

                labels[n] = genome_n
                node_id_to_group_idx[n] = n

                n += 1

    # ok, now we have a pile of k-mer vectors of size 4**args.ksize;
    # output in numpy format.

    # first, make a consistently ordered list of all k-mers, and convert
    # them into hashes.
    all_kmers = make_all(args.ksize)
    all_kmer_hashes = list(set([hash_murmur(i) for i in all_kmers]))
    all_kmer_hashes.sort()

    # now, build a matrix of GROUP_N rows x 4**ksize columns, where each
    # row will be the set of k-mer abundances associated with each group.
    V = numpy.zeros((len(group_info), 4 ** args.ksize), dtype=numpy.uint16)
    for i, n in enumerate(group_info):
        mh = group_info[n]
        vec = dict(mh.get_mins(with_abundance=True))
        vec = [vec.get(hashval, 0) for hashval in all_kmer_hashes]
        vec = numpy.array(vec)
        V[i] = vec

    # save!
    print("saving matrix of size {} to {}".format(str(V.shape), args.output))
    with open(args.output, "wb") as fp:
        numpy.save(fp, V)

    with open(args.output + ".labels", "wb") as fp:
        dump(labels, fp)

    with open(args.output + ".node_ids", "wb") as fp:
        pickle.dump(node_id_to_group_idx, fp)

    with open(args.output + ".node_mh", "wb") as fp:
        pickle.dump(group_ident, fp)


if __name__ == "__main__":
    main()
