#! /usr/bin/env python
import argparse
import os
import sys
import time
import pickle
from collections import defaultdict
import sqlite3

import screed
import sourmash

from ..utils.logging import notify, error
from . import MPHF_KmerIndex, hash_sequence
from .catlas import CAtlas
from . import search_utils


def main(argv):
    """\
    Query a catlas with a sequence (read, contig, or genome), and retrieve
    cDBG node IDs and MinHash signatures for the matching unitigs in the graph.
    """

    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument("cdbg_prefix", help="cdbg prefix")
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("output")
    p.add_argument("--query", help="query sequences", nargs="+")
    # @CTB add to config file for sgc
    p.add_argument(
        "-k", "--ksize", default=31, type=int, help="k-mer size (default: 31)"
    )
    p.add_argument('--query-is-dna', help='translate query to protein as well',
                   action='store_true')
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

    unitigs_db = os.path.join(args.cdbg_prefix, 'bcalm.unitigs.db')

    # get a single ksize for query
    prot_ksize = int(args.ksize)
    prot_mh = sourmash.MinHash(n=0, scaled=100, ksize=prot_ksize,
                               is_protein=True)

    record_kmers = defaultdict(set)
    all_kmers = set()
    
    if args.query_is_dna:
        add_as_protein = False
    else:
        add_as_protein = True

    for filename in args.query:
        print(f"Reading query from '{filename}'")

        screed_fp = screed.open(filename)
        for record in screed_fp:
            these_hashes = prot_mh.seq_to_hashes(record.sequence,
                                                 is_protein=add_as_protein)
            these_hashes = set(these_hashes)
            record_kmers[record.name] = these_hashes
            all_kmers.update(these_hashes)

    # iterate over unitigs next
    ## now unitigs... for all of the unitigs (which are DNA),
    ## first: translate the unitig into protein (up to six sequences),
    ## second: decompose into k-mers, save k-mers
    ## third: look for overlaps with query_kmers

    matching_cdbg = defaultdict(set)

    db = sqlite3.connect(unitigs_db)
    for n, record in enumerate(search_utils.contigs_iter_sqlite(db)):

        # translate into protein sequences
        seq = record.sequence

        unitig_hashes = prot_mh.seq_to_hashes(record.sequence)
        unitig_hashes = set(unitig_hashes)

        # do we have an overlap with any query??
        if unitig_hashes & all_kmers:
            # yes, match!
            cdbg_node = int(record.name)

            # iterate over all queries... slow.
            for query_name, query_hashes in record_kmers.items():
                # match to this record?
                if unitig_hashes & query_hashes:
                    matching_cdbg[query_name].add(cdbg_node)

        screed_fp.close()

    # ok, last iteration? expand neighborhoods.
    records_to_cdbg = {}
    cdbg_to_records = defaultdict(set)
    for query_name, cdbg_nodes in matching_cdbg.items():
        dominators = set()
        for cdbg_node in cdbg_nodes:
            dominators.add(catlas.cdbg_to_layer1[cdbg_node])

        print(f"got {len(dominators)} dominators for {record.name[:15]}")

        shadow = catlas.shadow(dominators)
        print(f"got {len(shadow)} cdbg_nodes under {len(dominators)} dominators")

        records_to_cdbg[('XXX', record.name)] = shadow
        for cdbg_node in shadow:
            cdbg_to_records[cdbg_node].add((filename, record.name))

    with open(outfile, "wb") as fp:
        print(f"saving pickled index to '{outfile}'")
        pickle.dump((args.catlas_prefix, records_to_cdbg, cdbg_to_records), fp)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))