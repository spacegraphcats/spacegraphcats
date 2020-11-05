#! /usr/bin/env python
import sys
import argparse
import time
from collections import defaultdict
import csv

import spacegraphcats
from spacegraphcats.utils.logging import notify, error, debug
from spacegraphcats.search import search_utils
from spacegraphcats.search import MPHF_KmerIndex, hash_sequence
from spacegraphcats.search.catlas import CAtlas
import screed


def main():
    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("query")
    p.add_argument("output")
    p.add_argument(
        "-k", "--ksize", default=31, type=int, help="k-mer size (default: 31)"
    )
    args = p.parse_args()

    # load catlas DAG
    catlas = CAtlas(args.catlas_prefix, load_sizefile=True)
    notify("loaded {} nodes from catlas {}", len(catlas), args.catlas_prefix)
    notify("loaded {} layer 1 catlas nodes", len(catlas.layer1_to_cdbg))

    ki_start = time.time()
    kmer_idx = MPHF_KmerIndex.from_catlas_directory(args.catlas_prefix)
    notify(
        "loaded {} k-mers in index ({:.1f}s)", len(kmer_idx), time.time() - ki_start,
    )

    # calculate the k-mer sizes for each catlas node.
    catlas.decorate_with_index_sizes(kmer_idx)

    notify(f"reading sequences from {args.query}")
    total = []
    for record in screed.open(args.query):
        total.append(record.sequence)

    seq = "".join(total)
    notify(f"...loaded {len(seq)} bp")

    notify("taking first 100k")
    seq = seq[:100000]

    hashes = hash_sequence(seq, args.ksize)
    notify(f"{len(hashes)} hashes total")

    # first, go through and get cDBG node IDs for each position
    cdbg_ids = []
    cdbg_id_kmers = defaultdict(set)

    n_none = 0
    n_total = 0
    for hash in hashes:
        n_total += 1

        # retrieve the cDBG ID for this k-mer
        cdbg_id = kmer_idx.get_cdbg_id(hash)
        cdbg_ids.append(cdbg_id)

        # track unique k-mers in this dataset that belong to this cDBG ID
        if cdbg_id is not None:
            cdbg_id_kmers[cdbg_id].add(hash)
        else:
            n_none += 1

    notify(f"done with first pass! {n_none} of {n_total} have no matches.")

    # now, for each layer1 dom node, get all cdbg nodes underneath, + sizes
    dom_sizes = {}

    # for each distinct cDBG ID,
    for cdbg_id in set(cdbg_ids):
        if cdbg_id is not None:
            # get the dom node
            dom_id = catlas.cdbg_to_layer1[cdbg_id]
            dom_cdbg_nodes = catlas.shadow([dom_id])
            assert dom_cdbg_nodes  # has to be at least 1!

            this_size = 0
            total_size = 0
            for cdbg_id in dom_cdbg_nodes:
                # number of k-mers from THIS query
                this_size += len(cdbg_id_kmers[cdbg_id])

                # number of k-mers from graph
                total_size += kmer_idx.get_cdbg_size(cdbg_id)

            dom_sizes[dom_id] = (this_size, total_size)

    # ok! do some counting -- for every k-mer, track dom node size, abund,
    # and containment of this data set in graph
    dom_node_sizes = []
    conts = []
    abunds = []
    for hash, cdbg_id in zip(hashes, cdbg_ids):
        if cdbg_id is not None:
            dom_id = catlas.cdbg_to_layer1[cdbg_id]
            dom_cdbg_nodes = catlas.shadow([dom_id])
            dom_node_sizes.append(len(dom_cdbg_nodes))

            # number of k-mers * unitig abund
            abund = catlas.weighted_kmer_sizes[cdbg_id]
            abunds.append(abund)

            query_size, domsize = dom_sizes[dom_id]
            containment = domsize / query_size
            if containment > 100:
                print("XYZ", query_size, domsize)
            conts.append(containment)
        else:
            dom_node_sizes.append(0)
            abunds.append(0)
            conts.append(0)

    with open(args.output, "wt") as fp:
        w = csv.writer(fp)
        w.writerow(["pos", "var", "neighbors", "abund"])
        for n, (cont, nbs, abund) in enumerate(zip(conts, dom_node_sizes, abunds)):
            w.writerow([n, cont, nbs, abund])

    return 0


if __name__ == "__main__":
    sys.exit(main())
