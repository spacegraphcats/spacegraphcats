#! /usr/bin/env python
"""
Look for catlas nodes that have few k-mers but many cDBG nodes underneath
them, as a likely sign of strain variation.
"""
import argparse
import os
import sys
import csv

import screed
from spacegraphcats.cdbg import hash_sequence
from .catlas import CAtlas
from . import search_utils


def main(args=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("query")
    p.add_argument("output")
    p.add_argument("--minsize", type=float, default=100)
    p.add_argument("--maxsize", type=float, default=10000)
    p.add_argument(
        "-k", "--ksize", default=31, type=int, help="k-mer size (default: 31)"
    )
    args = p.parse_args(args)

    print("minsize: {:g}".format(args.minsize))
    print("maxsize: {:g}".format(args.maxsize))

    basename = os.path.basename(args.catlas_prefix)
    # load catlas DAG
    catlas = CAtlas(args.catlas_prefix)
    top_node_id, dag, dag_levels = catlas.root, catlas.children, catlas.levels

    print("loaded {} nodes from catlas {}".format(len(dag), catlas))

    # load mapping between dom nodes and cDBG/graph nodes:
    layer1_to_cdbg = catlas.layer1_to_cdbg
    print("loaded {} layer 1 catlas nodes".format(len(layer1_to_cdbg)))

    # calculate the cDBG shadow sizes for each catlas node.
    print("decorating catlas with shadow size info.")
    catlas.decorate_with_shadow_sizes()
    node_shadow_sizes = catlas.shadow_sizes

    # ...and load cdbg node sizes
    print("loading contig size info")
    cdbg_kmer_sizes, cdbg_weighted_kmer_sizes = search_utils.load_cdbg_size_info(
        args.catlas_prefix
    )

    # decorate catlas with cdbg node sizes underneath them
    print("decorating catlas with contig size info.")
    (
        node_kmer_sizes,
        node_weighted_kmer_sizes,
    ) = search_utils.decorate_catlas_with_kmer_sizes(
        layer1_to_cdbg, dag, dag_levels, cdbg_kmer_sizes, cdbg_weighted_kmer_sizes
    )

    # load k-mer index, query, etc. etc.
    kmer_idx = search_utils.load_kmer_index(args.catlas_prefix)

    query_kmers = set()
    for record in screed.open(args.query):
        query_kmers.update(hash_sequence(record.sequence, args.ksize))

    print("got {}".format(len(query_kmers)))

    # construct dict cdbg_id -> # of query k-mers
    cdbg_match_counts = kmer_idx.count_cdbg_matches(query_kmers)

    total_match_kmers = sum(cdbg_match_counts.values())
    f_found = total_match_kmers / len(query_kmers)
    print("=> containment: {:.1f}%".format(f_found * 100))
    print("done loading & counting query k-mers in cDBG.")

    total_kmers_in_cdbg_matches = 0
    for cdbg_id in set(cdbg_match_counts.keys()):
        total_kmers_in_cdbg_matches += kmer_idx.get_cdbg_size(cdbg_id)

    cdbg_sim = total_match_kmers / total_kmers_in_cdbg_matches
    print("cdbg match node similarity: {:.1f}%".format(cdbg_sim * 100))
    cdbg_min_overhead = (
        total_kmers_in_cdbg_matches - total_match_kmers
    ) / total_match_kmers
    print("min cdbg overhead: {}".format(cdbg_min_overhead))

    # calculate the cDBG matching k-mers sizes for each catlas node.
    catlas_match_counts = kmer_idx.count_catlas_matches(
        cdbg_match_counts, dag, dag_levels, layer1_to_cdbg
    )

    ### ok, the real work: look at articulation of cDBG graph.

    # find highest nodes with kmer size less than given max_size
    def find_terminal_nodes(node_id, max_size):
        node_list = set()
        for sub_id in dag[node_id]:
            size = node_kmer_sizes[sub_id]

            if size < max_size:
                node_list.add(sub_id)
            else:
                children = find_terminal_nodes(sub_id, max_size)
                node_list.update(children)

        return node_list

    print("finding terminal nodes for {}.".format(args.maxsize))

    terminal = find_terminal_nodes(top_node_id, args.maxsize)
    print("...got {}".format(len(terminal)))
    terminal = {n for n in terminal if node_kmer_sizes[n] > args.minsize}
    print(
        "...down to {} between {} and {} in size.".format(
            len(terminal), args.minsize, args.maxsize
        )
    )

    # now, go through all nodes and print out characteristics
    with open(args.output, "wt") as fp:
        w = csv.writer(fp)

        w.writerow(
            ["node_id", "contained", "n_kmers", "n_weighted_kmers", "shadow_size"]
        )
        for n in terminal:
            f_contained = catlas_match_counts.get(n, 0) / node_kmer_sizes[n]
            w.writerow(
                [
                    str(n),
                    str(f_contained),
                    str(node_kmer_sizes[n]),
                    str(node_weighted_kmer_sizes[n]),
                    str(node_shadow_sizes[n]),
                ]
            )


if __name__ == "__main__":
    main()
