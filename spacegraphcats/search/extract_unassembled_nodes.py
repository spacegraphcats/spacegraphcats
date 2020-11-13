#! /usr/bin/env python
"""
Extract information about catlas nodes after splitting them into "assembled"
and "not".

EXPERIMENTAL.
"""
import argparse
import os
import sys
import csv
import sqlite3

import sourmash
import screed

from spacegraphcats.cdbg import hash_sequence
from . import search_utils
from .catlas import CAtlas


def main(args=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("query")
    p.add_argument("output")
    p.add_argument("--contigs-db", required=True)
    p.add_argument("--threshold", default=0.0, type=float)
    p.add_argument("--minsize", default=0, type=int)
    p.add_argument(
        "-k", "--ksize", default=31, type=int, help="k-mer size (default: 31)"
    )
    args = p.parse_args(args)

    print("threshold: {:.3f}".format(args.threshold))

    contigs_db = sqlite3.connect(args.contigs_db)

    basename = os.path.basename(args.catlas_prefix)
    catlas = CAtlas(args.catlas_prefix)

    # load catlas DAG
    print("loaded {} nodes from catlas {}".format(len(catlas), catlas))
    print("loaded {} layer 1 catlas nodes".format(len(catlas.layer1_to_cdbg)))

    # calculate the cDBG shadow sizes for each catlas node.
    print("decorating catlas with shadow size info.")
    catlas.decorate_with_shadow_sizes()

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
        catlas.layer1_to_cdbg,
        catlas.children,
        catlas.levels,
        cdbg_kmer_sizes,
        cdbg_weighted_kmer_sizes,
    )

    # load k-mer index, query, etc. etc.
    kmer_idx = search_utils.load_kmer_index(args.catlas_prefix)

    query_kmers = set()
    for record in screed.open(args.query):
        query_kmers.update(hash_sequence(record.sequence, args.ksize))

    print("got {} k-mers from {}".format(len(query_kmers), args.query))

    # construct dict cdbg_id -> # of query k-mers
    cdbg_match_counts = kmer_idx.count_cdbg_matches(query_kmers)

    total_match_kmers = sum(cdbg_match_counts.values())
    f_found = total_match_kmers / len(query_kmers)
    print("=> containment: {:.1f}%".format(f_found * 100))
    print("done loading & counting query k-mers in cDBG.")

    if total_match_kmers == 0:
        print("no match k-mers!?")
        sys.exit(-1)

    # calculate the cDBG matching k-mers sizes for each catlas node.
    catlas_match_counts = kmer_idx.count_catlas_matches(cdbg_match_counts, catlas)

    ### ok, the real work: find nodes that have low # of k-mers in the query.
    def find_unassembled_nodes(node_id, threshold=0.0):
        node_list = set()
        for sub_id in catlas.children[node_id]:
            n_matched = catlas_match_counts.get(sub_id, 0)
            size = node_kmer_sizes[sub_id]

            f_assembled = n_matched / size

            # if the fraction of unassembled k-mers under this node is below
            # our threshold, KEEP the node. Otherwise, descend into children.
            if f_assembled <= threshold:
                node_list.add(sub_id)
            else:
                children = find_unassembled_nodes(sub_id, threshold)
                node_list.update(children)

        return node_list

    print("finding unassembled nodes for threshold {}.".format(args.threshold))

    terminal = find_unassembled_nodes(catlas.root, args.threshold)
    sum_kmers = sum([node_kmer_sizes[n] for n in terminal])
    sum_match_kmers = sum([catlas_match_counts.get(n, 0) for n in terminal])
    print("...got {} nodes, representing {} k-mers".format(len(terminal), sum_kmers))

    # now, go through all nodes and print out characteristics
    print("writing node info to {}".format(args.output + ".csv"))
    with open(args.output + ".csv", "wt") as fp:
        w = csv.writer(fp)

        w.writerow(
            [
                "node_id",
                "contained",
                "n_kmers",
                "n_weighted_kmers",
                "average_weight",
                "shadow_size",
            ]
        )
        for n in terminal:
            f_contained = catlas_match_counts.get(n, 0) / node_kmer_sizes[n]
            w.writerow(
                [
                    n,
                    "{:.3f}".format(f_contained),
                    node_kmer_sizes[n],
                    "{:.1f}".format(node_weighted_kmer_sizes[n]),
                    "{:.2f}".format(node_weighted_kmer_sizes[n] / node_kmer_sizes[n]),
                    catlas.shadow_sizes[n],
                ]
            )

    if args.minsize:
        print("minsize set: {}. filtering.".format(args.minsize))
        new_terminal = set()
        for n in terminal:
            if node_kmer_sizes[n] >= args.minsize:
                new_terminal.add(n)

        print(
            "removed {} nodes => {}".format(
                len(terminal) - len(new_terminal), len(new_terminal)
            )
        )
        terminal = new_terminal

    # build cDBG shadow ID list, tagged by parent catlas node.
    cdbg_id_to_node = {}
    for n in terminal:
        shadow = catlas.shadow([n])
        for cdbg_id in shadow:
            cdbg_id_to_node[cdbg_id] = n

    #### extract contigs
    print("extracting contigs & building a sourmash signature")

    # track results as signature
    contigs_mh = sourmash.MinHash(n=0, ksize=args.ksize, scaled=1000)

    total_bp = 0
    total_seqs = 0

    print("writing contigs to {}".format(args.output + ".fa"))
    outfp = open(args.output + ".fa", "wt")
    for n, record in enumerate(search_utils.contigs_iter_sqlite(contigs_db)):
        if n and n % 10000 == 0:
            offset_f = total_seqs / len(cdbg_id_to_node)
            print(
                "...at n {} ({:.1f}% of shadow)".format(total_seqs, offset_f * 100),
                end="\r",
            )

        # contig names == cDBG IDs
        contig_id = int(record.name)
        catlas_parent = cdbg_id_to_node.get(contig_id)
        if catlas_parent is None:
            continue

        outfp.write(">{} {}\n{}\n".format(record.name, catlas_parent, record.sequence))
        contigs_mh.add_sequence(record.sequence)

        # track retrieved sequences in a minhash
        total_bp += len(record.sequence)
        total_seqs += 1

    # done - got all contigs!
    print("")
    print("fetched {} contigs, {} bp.".format(total_seqs, total_bp))

    print("writing sig to {}".format(args.output + ".sig"))
    with open(args.output + ".sig", "wt") as fp:
        ss = sourmash.SourmashSignature(contigs_mh)
        sourmash.save_signatures([ss], fp)

    return 0


if __name__ == "__main__":
    sys.exit(main())
