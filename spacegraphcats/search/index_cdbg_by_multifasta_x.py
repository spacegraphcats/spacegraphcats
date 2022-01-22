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
from spacegraphcats.utils.counter_gather import CounterGather
from .catlas import CAtlas
from . import search_utils


MAX_SQLITE_INT = 2 ** 63 - 1
sqlite3.register_adapter(
    int, lambda x: hex(x) if x > MAX_SQLITE_INT else x)
sqlite3.register_converter(
    'integer', lambda b: int(b, 16 if b[:2] == b'0x' else 10))


def load_all_signature_metadata(db):
    c = db.cursor()
    c.execute("SELECT id, name, num, scaled, ksize, filename, is_dna, is_protein, is_dayhoff, is_hp, track_abundance, seed FROM sketches")
    for (id, name, num, scaled, ksize, filename, is_dna, is_protein, is_dayhoff, is_hp, track_abundance, seed) in c:
        yield id, name, filename, ksize, scaled, is_protein


def load_sketch(db, sketch_id):
    c2 = db.cursor()

    c2.execute("SELECT name, num, scaled, ksize, filename, is_dna, is_protein, is_dayhoff, is_hp, track_abundance, seed FROM sketches WHERE id=?", (sketch_id,))

    name, num, scaled, ksize, filename, is_dna, is_protein, is_dayhoff, is_hp, track_abundance, seed = c2.fetchone()

    mh = sourmash.MinHash(n=num, ksize=ksize, scaled=scaled, seed=seed, is_protein=is_protein, dayhoff=is_dayhoff, hp=is_hp, track_abundance=track_abundance)

    c2.execute("SELECT hashval FROM hashes WHERE sketch_id=?", (sketch_id,))

    for hashval, in c2:
        print('XXX', hashval)
        mh.add_hash(hashval)

    ss = sourmash.SourmashSignature(mh, name=name, filename=filename)
    return ss


def get_matching_name_filename(query_cursor, unitig_mh):
    query_cursor.execute("DROP TABLE IF EXISTS hash_query")
    query_cursor.execute("CREATE TEMPORARY TABLE hash_query (hashval INTEGER PRIMARY KEY)")
    for hashval in unitig_mh.hashes:
        query_cursor.execute("INSERT INTO hash_query (hashval) VALUES (?)", (hashval,))

    #overlap = False
    #query_cursor.execute("SELECT EXISTS (SELECT 1 FROM hashes WHERE hashes.hashval IN (SELECT hashval FROM hash_query))")
    #overlap, = query_cursor.fetchone()

    # do we have an overlap with any query at all??
    query_cursor.execute("SELECT DISTINCT sketches.name,sketches.filename FROM sketches,hashes WHERE sketches.id=hashes.sketch_id AND hashes.hashval IN (SELECT hashval FROM hash_query)")

    for name, filename in query_cursor:
        yield name, filename


def get_matching_hashes(query_cursor, unitig_mh):
    query_cursor.execute("DROP TABLE IF EXISTS hash_query")
    query_cursor.execute("CREATE TEMPORARY TABLE hash_query (hashval INTEGER PRIMARY KEY)")
    for hashval in unitig_mh.hashes:
        query_cursor.execute("INSERT INTO hash_query (hashval) VALUES (?)", (hashval,))

    query_cursor.execute("SELECT DISTINCT hashes.hashval FROM hashes,hash_query WHERE hashes.hashval=hash_query.hashval")

    for hashval in query_cursor:
        yield hashval


def get_hashes_by_sketch_name_filename(query_cursor, name, filename):
    query_cursor.execute("SELECT DISTINCT hashes.hashval FROM hashes,sketches where sketches.name=? and sketches.filename=?", (name, filename))

    for hashval in query_cursor:
        yield hashval


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

    unitigs_db = os.path.join(args.cdbg_prefix, 'bcalm.unitigs.db')

    ### done loading! now let's do the thing.

    # track all hashes by record, as well as file origin of query
    record_hashes = defaultdict(set)
    query_idx_to_name = {}
    name_to_query_idx = {}
    hashval_to_queries = defaultdict(set)

    # CTB: as a future optimization, we could use an LCA database here
    # (on disk, etc.)

    # read all the queries into memory.
    this_query_idx = 0
    query_idx_to_hashes = {}
    all_kmers = set()
    assert len(args.query) == 1
    first_sketch = None
    for filename in args.query:
        notify(f"Reading query from database '{filename}'")

        query_db = sqlite3.connect(filename,
                                   detect_types=sqlite3.PARSE_DECLTYPES)


        for sketch_id, name, filename, ksize, scaled, is_protein in load_all_signature_metadata(query_db):
            first_sketch = load_sketch(query_db, sketch_id)
            assert first_sketch.minhash.is_protein
            assert first_sketch.minhash.ksize == args.ksize
            prot_mh = first_sketch.minhash.copy_and_clear()
            #print('XXX', prot_mh.scaled, prot_mh.is_protein, prot_mh.ksize)
            break

        for sketch_id, name, filename, ksize, scaled, is_protein in load_all_signature_metadata(query_db):
            assert prot_mh.scaled == scaled
            assert prot_mh.ksize == ksize
            assert is_protein

            query_idx_to_name[this_query_idx] = (filename, name)
            name_to_query_idx[(filename, name)] = this_query_idx
            this_query_idx += 1

    # iterate over unitigs next, to assign possible matches to each
    # (this is the prefetch stage.)
    #
    ## for all of the unitigs (which are DNA),
    ## first: translate the unitig into protein (up to six sequences),
    ## second: decompose into k-mers, save k-mers
    ## third: look for overlaps with query_kmers, save info on matches.

    query_matches_by_cdbg = {}

    notify(f"Iterating through unitigs in '{unitigs_db}'")
    db = sqlite3.connect(unitigs_db)
    query_cursor = query_db.cursor()
    for n, record in enumerate(search_utils.contigs_iter_sqlite(db)):
        if n % 100 == 0:
            print('...', n, len(query_matches_by_cdbg))
            if n >= 1000:
                break
        # translate into protein sequences
        unitig_mh = prot_mh.copy_and_clear()
        unitig_mh.add_sequence(record.sequence, force=True)

        if not unitig_mh:         # empty? skipme.
            continue

        matching_query_idx = set()
        for (name, filename) in get_matching_name_filename(query_cursor, unitig_mh):
            query_idx = name_to_query_idx[(filename, name)]
            matching_query_idx.add(query_idx)

        if matching_query_idx:
            # save match!
            cdbg_id = int(record.name)
            query_matches_by_cdbg[cdbg_id] = matching_query_idx

    # now, do some kind of summarization.
    #
    # options:
    # (1) any match whatsoever, promoted to neighborhood (OG mode)
    # (2) gather-filtered matches by cDBG node, promoted to neighborhood
    # (3) gather-filtered matches by neighborhood
    # (4) gather-filtered matches by entire graph - NOT IMPLEMENTED

    # cdbg_id => set([query_idx])
    filtered_query_matches_by_cdbg = {}

    notify(f"Doing filtering/promotion of matches: mode is {args.mode}")

    # (1) any match whatsoever, promoted to neighborhood
    if args.mode == "search+nbhd":
        # expand matches to neighborhood with no filtering
        for dom_id, shadow in catlas.layer1_to_cdbg.items():
            # collect query_idx from across neighborhood into a single set...
            this_query_idx = set()
            empty = set()
            for cdbg_id in shadow:
                mm = query_matches_by_cdbg.get(cdbg_id, empty)
                this_query_idx.update(mm)

            # ...and assign that set back to each cdbg ID
            if this_query_idx:
                for cdbg_id in shadow:
                    assert cdbg_id not in filtered_query_matches_by_cdbg
                    filtered_query_matches_by_cdbg[cdbg_id] = this_query_idx

    # (2) gather-filtered matches by cDBG node, promoted to neighborhood
    elif args.mode == "gather+cdbg":
        # look only at cDBG nodes with matches.
        for cdbg_id, this_query_idx in query_matches_by_cdbg.items():
            # for each node,

            # retrieve sequence
            record = search_utils.get_contigs_by_cdbg_sqlite(db, [cdbg_id])
            record = list(record)[0]

            # convert to hashes
            unitig_mh = prot_mh.copy_and_clear()
            unitig_mh.add_sequence(record.sequence, force=True)

            # get all matching kmers
            matching_kmers = set(get_matching_hashes(query_cursor, unitig_mh))
            assert matching_kmers

            # create gather counter with query == unitig hashes that match
            counter = CounterGather(matching_kmers)

            # retrieve matches by query_idx and fill counter
            for query_idx in this_query_idx:
                filename, name = query_idx_to_name[query_idx]
                hashes = set(get_hashes_by_sketch_name_filename(query_cursor, name, filename))
                counter.add(hashes, query_idx)

            # filter matches by min-set-cov
            this_filtered_query_idx = counter.do_full_gather()

            # update query_matches_by_cdbg with filtered indices
            query_matches_by_cdbg[cdbg_id] = set(this_filtered_query_idx)

        # now, expand matches to neighborhood
        for dom_id, shadow in catlas.layer1_to_cdbg.items():
            this_query_idx = set()

            # collect query_idx into a single set...
            empty = set()
            for cdbg_id in shadow:
                mm = query_matches_by_cdbg.get(cdbg_id, empty)
                this_query_idx.update(mm)

            if this_query_idx:
                # ...and assign that set back to each cdbg ID
                for cdbg_id in shadow:
                    assert cdbg_id not in filtered_query_matches_by_cdbg
                    filtered_query_matches_by_cdbg[cdbg_id] = this_query_idx

    # (3) gather-filtered matches by neighborhood
    elif args.mode == "gather+nbhd":
        # collect cdbg_ids by dominator; this involves iterating across
        # all dominators, which could be optimized...

        for dom_id, shadow in catlas.layer1_to_cdbg.items():
            # collect matching query_idx for all cDBG IDs under this dominator
            dom_matching_query_idx = set()
            for cdbg_id in shadow:
                if cdbg_id in query_matches_by_cdbg:
                    m = query_matches_by_cdbg[cdbg_id]
                    dom_matching_query_idx.update(m)

            # skip dominators that have no cDBG nodes with matches
            if not dom_matching_query_idx:
                continue

            # now, calculate & aggregate unitig hashes across dominator
            dom_mh = prot_mh.copy_and_clear()
            for record in search_utils.get_contigs_by_cdbg_sqlite(db, shadow):
                dom_mh.add_sequence(record.sequence, force=True)

            # get all matching kmers
            matching_kmers = set(get_matching_hashes(query_cursor, dom_mh))
            assert matching_kmers

            # create gather counter with query == dom hashes that match
            counter = CounterGather(matching_kmers)

            # build set of matches -
            for query_idx in dom_matching_query_idx:
                filename, name = query_idx_to_name[query_idx]
                hashes = set(get_hashes_by_sketch_name_filename(query_cursor, name, filename))
                counter.add(hashes, query_idx)

            # filter matches by min-set-cov
            this_filtered_query_idx = counter.do_full_gather()

            # store by cdbg_id
            for cdbg_id in shadow:
                # note: dominators should be disjoint, so we can directly
                # assign set rather than updating! but assert anyway.
                assert cdbg_id not in filtered_query_matches_by_cdbg
                filtered_query_matches_by_cdbg[cdbg_id] = \
                    set(this_filtered_query_idx)
    else:
        # (4) gather-filtered matches by entire graph NOT IMPLEMENTED yet ;).
        assert 0

    notify('...done!')

    # ok, last iteration: make output data structures by resolving
    # query_idx to (query_filename, query_name)

    records_to_cdbg = defaultdict(set)
    cdbg_to_records = {}
    for cdbg_id, query_idx_set in filtered_query_matches_by_cdbg.items():
        matches = set()
        for query_idx in query_idx_set:
            # @CTB SQL retrieve
            query_filename, query_name = query_idx_to_name[query_idx]

            matches.add((query_filename, query_name))
            records_to_cdbg[(query_filename, query_name)].add(cdbg_id)

        cdbg_to_records[cdbg_id] = matches

    if not records_to_cdbg:
        notify("WARNING: nothing in query matched to cDBG. Saving empty dictionaries.")

    # done! output.
    with open(outfile, "wb") as fp:
        notify(f"saving pickled index to '{outfile}'")
        pickle.dump((args.catlas_prefix, records_to_cdbg, cdbg_to_records), fp)
        notify(f"saved {len(records_to_cdbg)} query names with cDBG node mappings (of {len(query_idx_to_name)} queries total)")
        n_cdbg_match = len(cdbg_to_records)
        n_cdbg_total = len(catlas.cdbg_to_layer1)
        notify(f"saved {n_cdbg_match} cDBG IDs (of {n_cdbg_total} total; {n_cdbg_match / n_cdbg_total * 100:.1f}%) with at least one query match")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
