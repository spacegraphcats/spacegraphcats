#! /usr/bin/env python
"""
Choose catlas subtrees between no larger than --maxsize and no smaller
than --minsize in number of k-mers, and extract summary information on
abundances of kmers within the subtrees.
"""
import argparse
import os
import sys
import numpy
import pickle
import sqlite3

import sourmash
from sourmash.minhash import hash_murmur

from .catlas import CAtlas
from . import search_utils


def make_all(ksize):
    DNA = "ACGT"

    x = []

    def add(sofar, n):
        if n == 0:
            x.append("".join(sofar))
        else:
            for ch in DNA:
                add(sofar + ch, n - 1)

    add("", ksize)
    return x


def partition_catlas(catlas, max_size):
    roots = []

    def partition_recursive(node):
        if catlas.shadow_sizes[node] > max_size and len(catlas.children[node]) > 0:
            for u in catlas.children[node]:
                partition_recursive(u)
        else:
            roots.append(node)

    partition_recursive(catlas.root)
    return roots


def compute_matrix(group_info, group_ident, ksize, output):
    # first, make a consistently ordered list of all k-mers, and convert
    # them into hashes.
    all_kmers = make_all(ksize)
    all_kmer_hashes = list(set([hash_murmur(i) for i in all_kmers]))
    all_kmer_hashes.sort()

    # now, build a matrix of GROUP_N rows x 4**ksize columns, where each
    # row will be the set of k-mer abundances associated with each group.
    print("creating", len(group_info), 4 ** ksize)
    V = numpy.zeros((len(group_info), 4 ** ksize), dtype=numpy.uint16)
    node_id_to_group_idx = {}
    for i, n in enumerate(group_info):
        if i % 1000 == 0:
            print("...", i, len(group_info))
        mh = group_info[n]
        vec = mh.hashes
        vec = [vec.get(hashval, 0) for hashval in all_kmer_hashes]
        vec = numpy.array(vec)
        V[i] = vec

        node_id_to_group_idx[n] = i

    # save!
    print("saving matrix of size {} to {}".format(str(V.shape), output))
    with open(output, "wb") as fp:
        numpy.save(fp, V)

    with open(output + ".node_ids", "wb") as fp:
        pickle.dump(node_id_to_group_idx, fp)

    with open(output + ".node_mh", "wb") as fp:
        pickle.dump(group_ident, fp)


def main(args=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("output")
    p.add_argument("--contigs-db", required=True)
    p.add_argument("--maxsize", type=float, default=20000)
    p.add_argument("--minsize", type=float, default=5000)
    p.add_argument("--min-abund", type=float, default=0)
    p.add_argument("-k", "--ksize", default=5, type=int, help="k-mer size for vectors")
    p.add_argument("--scaled", type=int, default=1000)
    args = p.parse_args(args)

    print("minsize: {:g}".format(args.minsize))
    print("maxsize: {:g}".format(args.maxsize))
    print("ksize: {}".format(args.ksize))
    print("min_abund: {}".format(args.min_abund))

    catlas = CAtlas(args.catlas_prefix, load_sizefile=True, min_abund=args.min_abund)
    catlas.decorate_with_shadow_sizes()

    # everything is loaded!

    # find highest nodes with kmer size less than given max_size
    print("finding terminal nodes for {}.".format(args.maxsize))
    nodes = partition_catlas(catlas, args.maxsize)

    nodes = {n for n in nodes if catlas.kmer_sizes[n] > args.minsize}

    print(
        "{} nodes between {} and {} in k-mer size".format(
            len(nodes), args.minsize, args.maxsize
        )
    )
    print(
        "containing {} level1 nodes of {} total".format(
            len(catlas.shadow(nodes)), sum(map(len, catlas.layer1_to_cdbg.values()))
        )
    )

    node_kmers = sum([catlas.kmer_sizes[n] for n in nodes])
    total_kmers = catlas.kmer_sizes[catlas.root]
    print(
        "containing {} kmers of {} total ({:.1f}%)".format(
            node_kmers, total_kmers, node_kmers / total_kmers * 100
        )
    )

    # now build cdbg -> subtree/group ID

    cdbg_to_group = {}
    for n in nodes:
        shadow = catlas.shadow([n])
        for cdbg_id in shadow:
            # TODO remove cdbg vertices with no kmers
            # for cdbg_id in catlas.layer1_to_cdbg[level1_node]:
            # if cdbg_id in catlas.kmer_sizes:
            assert cdbg_id not in cdbg_to_group
            cdbg_to_group[cdbg_id] = n

    # record group info - here we are using the MinHash class to track
    # k-mer abundances in group_info, as well as using group_ident to
    # to track k=31 MinHashes for identification of each group.
    group_info = {}
    group_ident = {}
    for n in nodes:
        group_info[n] = sourmash.MinHash(
            n=0, ksize=args.ksize, scaled=1, track_abundance=1
        )
        group_ident[n] = sourmash.MinHash(n=0, ksize=31, scaled=args.scaled)

    # aaaaaand iterate over contigs, collecting abundances from all contigs
    # in a group.
    contigs_db = sqlite3.connect(args.contigs_db)
    for record_n, record in enumerate(search_utils.contigs_iter_sqlite(contigs_db)):
        if record_n % 10000 == 0:
            print("...", record_n, end="\r")
        cdbg_id = int(record.name)
        group_id = cdbg_to_group.get(cdbg_id)

        # if this is under a node that meets minsize criteria, track:
        if group_id is not None:
            # keep/measure abundances! CTB are actually doing anything abund?
            mh = group_info[group_id]
            mh.add_sequence(record.sequence, True)

            # update group idents.
            group_ident[group_id].add_sequence(record.sequence, True)

    # ok, now we have a pile of k-mer vectors of size 4**args.ksize;
    # output in numpy format.
    compute_matrix(group_info, group_ident, args.ksize, args.output)

    return 0


if __name__ == "__main__":
    sys.exit(main())
