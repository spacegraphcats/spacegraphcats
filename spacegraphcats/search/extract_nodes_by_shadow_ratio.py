#! /usr/bin/env python
"""
Look for catlas nodes that have few k-mers but many cDBG nodes underneath
them, as a likely sign of strain variation.

EXPERIMENTAL.
"""
import argparse
import os
import sys
import math
import sqlite3

import sourmash

import screed
from .catlas import CAtlas
from spacegraphcats.search import search_utils


def main(args=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("output")
    p.add_argument("--contigs-db", required=True)
    p.add_argument("--minsize", type=float, default=100)
    p.add_argument("--maxsize", type=float, default=10000)
    p.add_argument("--keep-fraction", type=float, default=0.1)
    p.add_argument(
        "-k", "--ksize", default=31, type=int, help="k-mer size (default: 31)"
    )
    args = p.parse_args(args)

    print("minsize: {:g}".format(args.minsize))
    print("maxsize: {:g}".format(args.maxsize))

    contigs_db = sqlite3.connect(args.contigs_db)

    # load catlas DAG
    catlas = CAtlas(args.catlas_prefix, load_sizefile=True)
    print("loaded {} nodes from catlas {}".format(len(catlas), catlas))
    print("loaded {} layer 1 catlas nodes".format(len(catlas.layer1_to_cdbg)))

    # calculate the cDBG shadow sizes for each catlas node.
    print("decorating catlas with shadow size info.")
    catlas.decorate_with_shadow_sizes()

    # ok, the real work: look at articulation of cDBG graph.

    # find highest nodes with kmer size less than given max_size
    def find_terminal_nodes(node_id, max_size):
        node_list = set()
        for sub_id in catlas.children[node_id]:
            # shadow size
            size = catlas.kmer_sizes[sub_id]

            if size < max_size:
                node_list.add(sub_id)
            else:
                children = find_terminal_nodes(sub_id, max_size)
                node_list.update(children)

        return node_list

    print("finding terminal nodes for {}.".format(args.maxsize))

    terminal = find_terminal_nodes(catlas.root, args.maxsize)
    print("...got {}".format(len(terminal)))
    terminal = {n for n in terminal if catlas.kmer_sizes[n] > args.minsize}
    print(
        "...down to {} between {} and {} in size.".format(
            len(terminal), args.minsize, args.maxsize
        )
    )

    # now, go through and calculate ratios
    x = []
    for node_id in terminal:
        # calculate: how many k-mers per cDBG node?
        kmer_size = catlas.kmer_sizes[node_id]
        shadow_size = catlas.shadow_sizes[node_id]

        ratio = math.log(kmer_size, 2) - math.log(shadow_size, 2)

        # track basic info
        x.append((ratio, node_id, shadow_size, kmer_size))

    print("terminal node stats for maxsize: {:g}".format(args.maxsize))
    print("n tnodes:", len(terminal))
    print("total k-mers:", catlas.kmer_sizes[catlas.root])

    x.sort(reverse=True)
    for (k, v, a, b) in x[:10]:
        print("ratio: {:.3f}".format(2 ** k), "/ shadow size:", a, "/ kmers:", b)
    print("... eliding {} nodes".format(len(x) - 20))
    for (k, v, a, b) in x[-10:]:
        print("ratio: {:.3f}".format(2 ** k), "/ shadow size:", a, "/ kmers:", b)

    # keep the last keep-fraction (default 10%) for examination
    keep_sum_kmer = args.keep_fraction * catlas.kmer_sizes[catlas.root]
    sofar = 0
    keep_terminal = set()
    for (k, v, a, b) in reversed(x):
        sofar += b
        if sofar > keep_sum_kmer:
            break
        keep_terminal.add(v)

    print(f"keeping last {sofar} k-mers worth of nodes for examination.")

    # build cDBG shadow ID list.
    cdbg_shadow = catlas.shadow(keep_terminal)

    # track results as signature
    contigs_mh = sourmash.MinHash(n=0, ksize=args.ksize, scaled=1000)

    total_bp = 0
    total_seqs = 0

    outfp = open(args.output, "wt")
    for record in search_utils.get_contigs_by_cdbg_sqlite(contigs_db, cdbg_shadow):
        outfp.write(">{}\n{}\n".format(record.name, record.sequence))
        contigs_mh.add_sequence(record.sequence)

        # track retrieved sequences in a minhash
        total_bp += len(record.sequence)
        total_seqs += 1

    # done - got all contigs!
    print("")
    print("fetched {} contigs, {} bp.".format(total_seqs, total_bp))

    print("wrote contigs to {}".format(args.output))
    with open(args.output + ".sig", "wt") as fp:
        ss = sourmash.SourmashSignature(contigs_mh)
        sourmash.save_signatures([ss], fp)

    return 0


if __name__ == "__main__":
    sys.exit(main())
