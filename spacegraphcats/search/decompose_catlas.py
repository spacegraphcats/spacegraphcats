#! /usr/bin/env python
"""
Choose catlas subtrees between no larger than --maxsize and no smaller
than --minsize in number of k-mers, and extract summary information on
abundances of kmers within the subtrees.
"""
import argparse
import os
import sys
import gzip

from .catlas import CAtlas


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


def main(args=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("output_directory")
    p.add_argument("--maxsize", type=float, default=20000)
    p.add_argument("--minsize", type=float, default=5000)
    p.add_argument("--min-abund", type=float, default=0)
    p.add_argument("--scaled", type=int, default=1000)
    args = p.parse_args(args)

    os.mkdir(args.output_directory)

    print("minsize: {:g}".format(args.minsize))
    print("maxsize: {:g}".format(args.maxsize))
    print("min_abund: {}".format(args.min_abund))

    catlas = CAtlas(args.catlas_prefix, load_sizefile=True, min_abund=args.min_abund)
    catlas.decorate_with_shadow_sizes()

    # everything is loaded!

    # find highest nodes with kmer size less than given max_size
    print("finding terminal nodes for {}.".format(args.maxsize))
    nodes = partition_catlas(catlas, args.maxsize)

    nodes = {n for n in nodes if catlas.kmer_sizes[n] > args.minsize}

    shadow = catlas.shadow(nodes)

    print(
        "{} nodes between {} and {} in k-mer size".format(
            len(nodes), args.minsize, args.maxsize
        )
    )
    print(
        "containing {} level1 nodes of {} total".format(
            len(shadow), sum(map(len, catlas.layer1_to_cdbg.values()))
        )
    )

    node_kmers = sum([catlas.kmer_sizes[n] for n in nodes])
    total_kmers = catlas.kmer_sizes[catlas.root]
    print(
        "containing {} kmers of {} total ({:.1f}%)".format(
            node_kmers, total_kmers, node_kmers / total_kmers * 100
        )
    )

    print("outputting cDBG node IDs to directory {}".format(args.output_directory))

    for n in nodes:
        with gzip.open("{}/{}.txt.gz".format(args.output_directory, n), "wt") as fp:
            fp.write("\n".join([str(x) for x in catlas.shadow([n])]))


if __name__ == "__main__":
    main()
