"""
Annotate/index the cDBG with FASTA records from one or more query files,
in protein space.

This can be used, for example, to connect cDBG nodes to gene names.

In brief, this script
- reads query sequences from FASTA files
- records hashes from query sequences in protein space
- iterates over the cDBG unitigs and finds cDBG nodes that match to each
  query. (This is potentially a quadratic step.)
- optionally does filtering and/or expansion of annotations using min-set-cov/
  gather and neighborhoods/dominators.
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
from collections import defaultdict
import sqlite3

import screed
import sourmash

from ..utils.logging import notify, error
from .catlas import CAtlas
from . import search_utils


def prefetch(query_mh, annot_databases):
    if not query_mh:
        return

    results = []
    query = sourmash.SourmashSignature(query_mh)
    for annot_db in annot_databases:
        for match in annot_db.prefetch(query, 0):
            results.append(match.signature)

    return results


def gather(query_mh, matches):
    counter = sourmash.index.CounterGather(query_mh)
    for m in matches:
        counter.add(m)

    # build the min set cov
    results = []
    while 1:
        result = counter.peek(query_mh, 0)
        if result:
            (sr, intersect_mh) = result
            counter.consume(intersect_mh)
            results.append(sr.signature)
        else:
            break

    return results


def batch_process_dominators(dom_to_seq, prot_mh, annot_dblist):
    all_mh = prot_mh.copy_and_clear()
    for dom_id, records in dom_to_seq.items():
        for record in records:
            all_mh.add_sequence(record.sequence, force=True)

    print('XXX', len(all_mh.hashes))

    all_matches = prefetch(all_mh, annot_dblist)
    return sourmash.index.LinearIndex(_signatures=all_matches)


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
    p.add_argument('--mode', type=str, default="search+nbhd")
    p.add_argument('--query-by-file', default=False, action="store_true")

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
        error("must specify at least one query file using --query.")
        sys.exit(-1)

    if args.mode not in ("search+nbhd", "gather+cdbg", "gather+nbhd"):
        error("Please specify a valid match mode with --mode.")
        error("Valid modes are: 'search+nbhd', 'gather+cdbg', 'gather+nbhd'")
        sys.exit(-1)
    notify(f"Using query mode: {args.mode}")

    # make sure all of the query sequences exist.
    for filename in args.query:
        if not os.path.exists(filename):
            error("query db file {} does not exist.", filename)
            sys.exit(-1)

    if args.query_by_file:
        notify("Aggregating queries by file because --query-by-file specified.")

    else:
        notify("Using individual records from query files.")

    # get a single ksize for query
    prot_ksize = int(args.ksize)
    notify(f"Using protein k-mer size {prot_ksize}")

    # translate query, or not?
    if args.query_is_dna:
        assert False # @CTB
        add_as_protein = False
        notify("Translating queries from DNA into protein.")
    else:
        add_as_protein = True
        notify("Queries are proteins; no translation occurring.")

    # load catlas DAG
    catlas = CAtlas(args.cdbg_prefix, args.catlas_prefix)
    notify(f"loaded {len(catlas)} nodes from catlas {args.catlas_prefix}")
    notify(f"loaded {len(catlas.layer1_to_cdbg)} layer 1 catlas nodes")

    unitigs_db_path = os.path.join(args.cdbg_prefix, 'bcalm.unitigs.db')
    unitigs_db = sqlite3.connect(unitigs_db_path)

    ### done loading! now let's do the thing.

    # track all hashes by record, as well as file origin of query
    record_hashes = defaultdict(set)

    # CTB: as a future optimization, we could use an LCA database here
    # (on disk, etc.)

    assert args.mode in ("search+nbhd", "gather+nbhd")

    # read all the queries into memory.
    this_query_idx = 0
    query_idx_to_hashes = {}
    all_kmers = set()
    first_sketch = None

    assert len(args.query) == 1
    annot_database = args.query[0]
    notify(f"Reading annotation database: '{annot_database}'")
    annot_db = sourmash.load_file_as_index(annot_database)
    annot_dblist = [annot_db]

    # get protein MinHash object to use. @CTB use select etc etc here
    for ss in annot_dblist[0].signatures():
        prot_mh = ss.minhash.copy_and_clear()
        assert prot_mh.is_protein
        assert prot_mh.ksize == args.ksize
        break

    cdbg_to_records = {}

    # collect cdbg_ids by dominator; this involves iterating across
    # all dominators, which could maybe be optimized...
    # then aggregate into one set of hashes, and run gather.

    print(f'Iterating over {len(catlas.layer1_to_cdbg)} dominators.')
    batch = {}
    for n, (dom_id, shadow) in enumerate(catlas.layer1_to_cdbg.items()):
        if n % 1000 == 0:
            print('...', n)

        if (n % 1000 == 0 or n == len(catlas.layer1_to_cdbg) - 1) and batch:
            matches_db = batch_process_dominators(batch, prot_mh,
                                                  annot_dblist)

            if matches_db:
                # now, process each dom_id individually.
                for dom_id in batch:
                    shadow = catlas.layer1_to_cdbg[dom_id]
                    records = batch[dom_id]

                    dom_mh = prot_mh.copy_and_clear()
                    for record in records:
                        dom_mh.add_sequence(record.sequence, force=True)

                    # get overlapping signatures
                    matches = prefetch(dom_mh, [matches_db])
                    if not matches:
                        continue

                    matching_annots = set()

                    # filter? do a gather, if requested
                    if args.mode == 'gather+nbhd':
                        results = gather(dom_mh, matches)

                        for m in results:
                            matching_annots.add((m.filename, m.name))
                    else:
                        # no filter.
                        assert args.mode == 'search+nbhd'
                        for m in matches:
                            matching_annots.add((m.filename, m.name))

                    # store by cdbg_id
                    for cdbg_id in shadow:
                        # note: dominators should be disjoint, so we can directly
                        # assign set rather than updating! but assert anyway.
                        assert cdbg_id not in cdbg_to_records
                        cdbg_to_records[cdbg_id] = set(matching_annots)

            batch = {}
            shadows = {}

        records = []
        for record in search_utils.get_contigs_by_cdbg_sqlite(unitigs_db,
                                                              shadow):
            records.append(record)
        batch[dom_id] = records

    ###

    notify('...done!')

    # ok, last iteration: make output data structures by resolving
    # query_idx to (query_filename, query_name)

    records_to_cdbg = defaultdict(set)
    for cdbg_id, matching_annots in cdbg_to_records.items():
        for (filename, name) in matching_annots:
            records_to_cdbg[(filename, name)].add(cdbg_id)

    if not records_to_cdbg:
        notify("WARNING: nothing in query matched to cDBG. Saving empty dictionaries.")

    # done! output.
    with open(outfile, "wb") as fp:
        notify(f"saving pickled index to '{outfile}'")
        pickle.dump((args.catlas_prefix, records_to_cdbg, cdbg_to_records), fp)
        notify(f"saved {len(records_to_cdbg)} query names with cDBG node mappings")
        n_cdbg_match = len(cdbg_to_records)
        n_cdbg_total = len(catlas.cdbg_to_layer1)
        notify(f"saved {n_cdbg_match} cDBG IDs (of {n_cdbg_total} total; {n_cdbg_match / n_cdbg_total * 100:.1f}%) with at least one query match")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
