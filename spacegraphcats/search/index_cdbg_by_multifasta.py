#! /usr/bin/env python
"""
Produce an index that connects cDBG nodes to matching records in a FASTA file.

This is primarily used for annotating cDBG neighborhoods with FASTA records.

In brief,
* load catlas & k-mer index into cDBG
* for every record in one or more query FASTA files,
  - find all cDBG nodes that match to the k-mers in that record
  - find all dominators that contain those cDBG nodes, and expand cDBG IDs
    to all nodes under those dominators
  - build dictionary 'cdbg_to_records',
    'cdbg_id' => set( (query_file, record_name) )
  - build dictionary 'records_to_cdbg',
    ('query_file, record_name') => set( cdbg_ids )

* save constructed dictionaries to a pickle file as index.

See index_cdbg_by_multifasta_x for a protein-space version of this script.
"""
import argparse
import os
import sys
import time
import pickle
from collections import defaultdict

import screed

from ..utils.logging import notify, error
from . import MPHF_KmerIndex, hash_sequence
from .catlas import CAtlas


def main(argv):
    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument("cdbg_prefix", help="cdbg prefix")
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("output")
    p.add_argument("--query", help="query sequences", nargs="+")
    p.add_argument("-v", "--verbose", action="store_true")

    args = p.parse_args(argv)
    outfile = args.output

    if not args.query:
        print("must specify at least one query file using --query.")
        sys.exit(-1)

    # make sure all of the query sequences exist.
    for filename in args.query:
        if not os.path.exists(filename):
            error("query seq file {} does not exist.", filename)
            sys.exit(-1)

    # load catlas DAG
    catlas = CAtlas(args.cdbg_prefix, args.catlas_prefix)
    notify("loaded {} nodes from catlas {}", len(catlas), args.catlas_prefix)
    notify("loaded {} layer 1 catlas nodes", len(catlas.layer1_to_cdbg))

    # ...and kmer index.
    ki_start = time.time()
    kmer_idx = MPHF_KmerIndex.from_directory(args.cdbg_prefix)

    ksize = kmer_idx.ksize
    notify(f"Using ksize {ksize} from k-mer index.")
    notify("loaded {} k-mers in index ({:.1f}s)", len(kmer_idx), time.time() - ki_start)

    # calculate the k-mer sizes for each catlas node.
    catlas.decorate_with_index_sizes(kmer_idx)

    # use the same ksize as the kmer index.

    records_to_cdbg = {}
    n_records_found = 0
    cdbg_to_records = defaultdict(set)

    for filename in args.query:
        print(f"Reading from '{filename}'")
        screed_fp = screed.open(filename)
        for record in screed_fp:
            if len(record.sequence) < int(ksize):
                continue

            kmers = hash_sequence(record.sequence, ksize)
            cdbg_match_counts = kmer_idx.count_cdbg_matches(kmers)

            print(
                f"got {len(cdbg_match_counts)} cdbg nodes for {record.name[:15]} ({len(kmers)} kmers)"
            )

            dominators = set()
            for cdbg_node in cdbg_match_counts:
                dominators.add(catlas.cdbg_to_layer1[cdbg_node])

            print(f"got {len(dominators)} dominators for {record.name[:15]}")

            shadow = catlas.shadow(dominators)
            print(f"got {len(shadow)} cdbg_nodes under {len(dominators)} dominators")

            records_to_cdbg[(filename, record.name)] = shadow
            if shadow:
                n_records_found += 1
                for cdbg_node in shadow:
                    cdbg_to_records[cdbg_node].add((filename, record.name))

        screed_fp.close()

    if not records_to_cdbg:
        print(
            "WARNING: nothing in query matched to cDBG. Saving empty dictionaries.",
            file=sys.stderr,
        )

    with open(outfile, "wb") as fp:
        print(f"saving pickled index to '{outfile}'")
        pickle.dump((args.catlas_prefix, records_to_cdbg, cdbg_to_records), fp)
        print(
            f"saved {n_records_found} query names with cDBG node mappings (of {len(records_to_cdbg)} queries total)"
        )
        n_cdbg_match = len(cdbg_to_records)
        n_cdbg_total = len(catlas.cdbg_to_layer1)
        print(
            f"saved {n_cdbg_match} cDBG IDs (of {n_cdbg_total} total; {n_cdbg_match / n_cdbg_total * 100:.1f}%) with at least one query match"
        )

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
