#! /usr/bin/env python
"""
Annotate/index the cDBG with FASTA records from one or more query files,
in protein space.

This can be used, for example, to connect cDBG nodes to gene names.

In brief, this script
- reads query sequences from FASTA files
- records hashes from query sequences in protein space
- iterates over the cDBG unitigs and finds cDBG nodes that match to each
  query. (This is potentially a quadratic step.)
- produces mappings from query to cDBG nodes, and vice versa.
- saves to pickle file.

See index_cdbg_by_multifasta for a pure DNA-space version of this script.
"""
import argparse
import os
import sys
import pickle
from collections import defaultdict
import sqlite3

import screed
import sourmash

from ..utils.logging import notify, error
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

    p.add_argument(
        "-k", "--ksize", default=10, type=int,
        help="protein k-mer size (default: 10)"
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

    # get a single ksize for query
    prot_ksize = int(args.ksize)
    prot_mh = sourmash.MinHash(n=0, scaled=100, ksize=prot_ksize,
                               is_protein=True)
    print(f"Using protein k-mer size {prot_ksize}")

    # translate query, or not?
    if args.query_is_dna:
        add_as_protein = False
        print("Translating queries from DNA into protein.")
    else:
        add_as_protein = True
        print("Queries are proteins; no translation occurring.")

    # load catlas DAG
    catlas = CAtlas(args.cdbg_prefix, args.catlas_prefix)
    notify("loaded {} nodes from catlas {}", len(catlas), args.catlas_prefix)
    notify("loaded {} layer 1 catlas nodes", len(catlas.layer1_to_cdbg))

    unitigs_db = os.path.join(args.cdbg_prefix, 'bcalm.unitigs.db')

    ### done loading! now let's do the thing.

    # track all hashes by record, as well as file origin of query
    record_hashes = defaultdict(set)
    query_name_to_filename = {}
    all_kmers = set()
    
    # read all the queries into memory.
    for filename in args.query:
        print(f"Reading query from '{filename}'")

        screed_fp = screed.open(filename)
        for record in screed_fp:
            these_hashes = prot_mh.seq_to_hashes(record.sequence,
                                                 is_protein=add_as_protein)
            these_hashes = set(these_hashes)
            record_hashes[record.name] = these_hashes
            all_kmers.update(these_hashes)

            if record.name in query_name_to_filename:
                print(f"ERROR: record '{record.name}' from '{filename}' is a duplicate!", file=sys.stderr)
                print(f"...also appears in '{query_name_to_filename[record.name]}'", file=sys.stderr)
                print("Failing because this violates our assumptions in this script.", file=sys.stderr)
                sys.exit(-1)

            query_name_to_filename[record.name] = filename

        screed_fp.close()

    # iterate over unitigs next
    ## now unitigs... for all of the unitigs (which are DNA),
    ## first: translate the unitig into protein (up to six sequences),
    ## second: decompose into k-mers, save k-mers
    ## third: look for overlaps with query_kmers

    matching_cdbg = defaultdict(set)

    db = sqlite3.connect(unitigs_db)
    for n, record in enumerate(search_utils.contigs_iter_sqlite(db)):
        # translate into protein sequences
        unitig_hashes = prot_mh.seq_to_hashes(record.sequence)
        unitig_hashes = set(unitig_hashes)

        # do we have an overlap with any query??
        if unitig_hashes & all_kmers:
            # yes, match!
            cdbg_node = int(record.name)

            # iterate over all queries... this is where things get potentially
            # slow, b/c this is quadratic!
            for query_name, query_hashes in record_hashes.items():
                # match to this record?
                if unitig_hashes & query_hashes:
                    matching_cdbg[query_name].add(cdbg_node)

        screed_fp.close()

    # ok, last iteration? expand neighborhoods.
    records_to_cdbg = {}
    cdbg_to_records = defaultdict(set)
    for query_name, cdbg_nodes in matching_cdbg.items():
        query_filename = query_name_to_filename[query_name]
        dominators = set()
        for cdbg_node in cdbg_nodes:
            dominators.add(catlas.cdbg_to_layer1[cdbg_node])

        print(f"got {len(dominators)} dominators for {query_name[:15]}")

        shadow = catlas.shadow(dominators)
        print(f"got {len(shadow)} cdbg_nodes under {len(dominators)} dominators")

        records_to_cdbg[(query_filename, query_name)] = shadow
        for cdbg_node in shadow:
            cdbg_to_records[cdbg_node].add((filename, record.name))

    if not records_to_cdbg:
        print("WARNING: nothing in query matched to cDBG. Saving empty dictionaries.", file=sys.stderr)

    with open(outfile, "wb") as fp:
        print(f"saving pickled index to '{outfile}'")
        pickle.dump((args.catlas_prefix, records_to_cdbg, cdbg_to_records), fp)

    # @CTB remove.
    import pprint
    pprint.pprint(records_to_cdbg)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
