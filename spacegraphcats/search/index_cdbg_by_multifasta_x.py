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
- expands cDBG unitigs into neighborhood using level 1 dominators
- produces mappings from query to cDBG nodes, and vice versa.
- saves to pickle file.

See index_cdbg_by_multifasta for a pure DNA-space version of this script.

NOTE: unlike index_cdbg_by_multifasta, query records with no match are NOT
saved.
"""
import argparse
import os
import sys
import pickle
from collections import defaultdict, Counter
import sqlite3

import screed
import sourmash

from ..utils.logging import notify, error
from .catlas import CAtlas
from . import search_utils


class CounterGather:
    """
    A refactoring of sourmash.index.CounterGather to use sets of hashes.

    The public interface is `peek(...)` and `consume(...)` only.
    """
    def __init__(self, query_set):
        # track query
        self.orig_query_set = set(query_set)

        # track matches in a list & their locations
        self.setlist = []
        self.locations = []

        # ...and overlaps with query
        self.counter = Counter()

        # cannot add matches once query has started.
        self.query_started = 0

    def add(self, hashes, ident, require_overlap=True):
        "Add this set of hashes in as a potential match."
        if self.query_started:
            raise ValueError("cannot add more hashes to counter after peek/consume")

        # upon insertion, count & track overlap with the specific query.
        overlap = self.orig_query_set & hashes
        if overlap:
            i = len(self.setlist)

            self.counter[i] = len(overlap)
            self.setlist.append(hashes)
            self.locations.append(ident)
        elif require_overlap:
            raise ValueError("no overlap between query and signature!?")

    def peek(self, cur_query_set):
        "Get next 'gather' result for this database, w/o changing counters."
        self.query_started = 1

        # empty? nothing to search.
        counter = self.counter
        if not counter:
            return []

        setlist = self.setlist
        assert setlist

        if not cur_query_set:             # empty query? quit.
            return []

        if cur_query_set & self.orig_query_set != cur_query_set:
            raise ValueError("current query not a subset of original query")

        # Find the best match -
        most_common = counter.most_common()
        dataset_id, match_size = most_common[0]

        # pull match and location.
        match_set = setlist[dataset_id]

        # calculate containment
        cont = len(cur_query_set & match_set) / len(match_set)
        assert cont

        # calculate intersection of this "best match" with query.
        intersect_set = cur_query_set & match_set
        location = self.locations[dataset_id]

        # build result & return intersection
        return cont, match_set, location, intersect_set

    def consume(self, intersect_set):
        "Remove the given hashes from this counter."
        self.query_started = 1

        if not intersect_set:
            return

        setlist = self.setlist
        counter = self.counter

        most_common = counter.most_common()

        # Prepare counter for finding the next match by decrementing
        # all hashes found in the current match in other datasets;
        # remove empty datasets from counter, too.
        for (dataset_id, _) in most_common:
            remaining_set = setlist[dataset_id]
            intersect_count = len(intersect_set & remaining_set)
            if intersect_count:
                counter[dataset_id] -= intersect_count
                if counter[dataset_id] == 0:
                    del counter[dataset_id]


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
    query_idx_to_name = {}
    hashval_to_queries = defaultdict(set)

    # @CTB: as a future optimization, we could use an LCA database here
    # (on disk, etc.)
    # read all the queries into memory.
    this_query_idx = 0
    all_kmers = set()
    for filename in args.query:
        print(f"Reading query from '{filename}'")

        screed_fp = screed.open(filename)
        for record in screed_fp:
            these_hashes = prot_mh.seq_to_hashes(record.sequence,
                                                 is_protein=add_as_protein)

            query_idx_to_name[this_query_idx] = (filename, record.name)
            for hashval in these_hashes:
                hashval_to_queries[hashval].add(this_query_idx)

            all_kmers.update(these_hashes)
            
            this_query_idx += 1

        screed_fp.close()

    # iterate over unitigs next
    ## now unitigs... for all of the unitigs (which are DNA),
    ## first: translate the unitig into protein (up to six sequences),
    ## second: decompose into k-mers, save k-mers
    ## third: look for overlaps with query_kmers

    matching_cdbg = defaultdict(set)

    print(f"Iterating through unitigs in '{unitigs_db}'")
    db = sqlite3.connect(unitigs_db)
    for n, record in enumerate(search_utils.contigs_iter_sqlite(db)):
        # translate into protein sequences
        unitig_hashes = prot_mh.seq_to_hashes(record.sequence)

        # do we have an overlap with any query??
        matching_kmers = set(unitig_hashes) & all_kmers
        matching_query_idx = set()
        for hashval in matching_kmers:
            matching_query_idx.update(hashval_to_queries[hashval])

        if matching_query_idx:
            # yes, match!
            cdbg_node = int(record.name)

            # @CTB: here we probably want to track hashes, instead of
            # nodes.
            for query_idx in matching_query_idx:
                matching_cdbg[query_idx].add(cdbg_node)

        screed_fp.close()

    print('...done!')

    # @CTB: we might want to (optionally) expand neighborhoods first,
    # to get list/set of hashes.

    print('Expanding neighborhoods:')

    # ok, last iteration? expand neighborhoods.
    records_to_cdbg = {}
    cdbg_to_records = defaultdict(set)
    for query_idx, cdbg_nodes in matching_cdbg.items():
        query_filename, query_name = query_idx_to_name[query_idx]
        dominators = set()
        for cdbg_node in cdbg_nodes:
            dominators.add(catlas.cdbg_to_layer1[cdbg_node])

        print(f"got {len(dominators)} dominators for {query_name[:15]}")

        shadow = catlas.shadow(dominators)
        print(f"got {len(shadow)} cdbg_nodes under {len(dominators)} dominators")

        records_to_cdbg[(query_filename, query_name)] = shadow

        # @CTB: here, we want to do something with hashes?
        for cdbg_node in shadow:
            cdbg_to_records[cdbg_node].add((filename, record.name))

    if not records_to_cdbg:
        print("WARNING: nothing in query matched to cDBG. Saving empty dictionaries.", file=sys.stderr)

    with open(outfile, "wb") as fp:
        print(f"saving pickled index to '{outfile}'")
        pickle.dump((args.catlas_prefix, records_to_cdbg, cdbg_to_records), fp)
        print(f"saved {len(records_to_cdbg)} query names with cDBG node mappings (of {this_query_idx + 1} queries total)")
        n_cdbg_match = len(cdbg_to_records)
        n_cdbg_total = len(catlas.cdbg_to_layer1)
        print(f"saved {n_cdbg_match} cDBG IDs (of {n_cdbg_total} total; {n_cdbg_match / n_cdbg_total * 100:.1f}%) with at least one query match")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
