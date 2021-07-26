#! /usr/bin/env python
"""
Sort the bcalm unitigs.fa output (a cDBG) into deterministic order.
Produces output consumed by `bcalm_to_gxt`.

Also saves all of the graph information for later perusal.

Also outputs a sourmash scaled=1000 signature for the input unitigs.
"""
import sys
import argparse
import pickle
import collections
import sqlite3

import sourmash
from spacegraphcats.search.search_utils import my_fasta_iter


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("gfa_unitigs", help="input GFA file for cDBG")
    parser.add_argument("sqlite_db_out",
                        help="output SQLite databsae for unitigs")
    parser.add_argument("mapping_pickle_out",
                        help="output pickle file for graph structure")
    parser.add_argument("--sig")
    parser.add_argument("-k", "--ksize", type=int, default=31,
                        help="k-mer size used to build cDBG")
    parser.add_argument("--seed", type=int, default=42,
                        help="minimizer seed for sorting contigs")
    args = parser.parse_args(argv)

    unitigs = args.gfa_unitigs
    ksize = args.ksize

    # this MinHash is used to find a minimizer for each contig, so that
    # can sort them deterministically.
    empty_mh = sourmash.MinHash(n=1, ksize=ksize, seed=args.seed)

    ###
    db = sqlite3.connect(args.sqlite_db_out)
    cursor = db.cursor()
    cursor.execute("PRAGMA synchronous='OFF'")
    cursor.execute("PRAGMA locking_mode=EXCLUSIVE")
    cursor.execute(
        "CREATE TABLE sequences (id INTEGER, sequence TEXT, offset INTEGER PRIMARY KEY, min_hashval TEXT, abund FLOAT)"
    )

    # read & find minimum hash:
    total_bp = 0
    neighbors = collections.defaultdict(set)

    # record input k-mers in a minhash
    in_mh = sourmash.MinHash(0, ksize, scaled=1000)

    print(f"loading GFA unitigs from {unitigs}...")

    unitigs_fp = open(unitigs, "rt")

    def iter_lines_and_pos(fp):
        pos = fp.tell()
        line = fp.readline()
        while line:
            yield line, pos
            pos = fp.tell()
            line = fp.readline()
    
    last_pos = 0
    k_overlap_str = f"{ksize-1}M"
    for n, (line, offset) in enumerate(iter_lines_and_pos(unitigs_fp)):
        # NOTE: here, 'offset' is used as a unique identifier for 'S'
        # records, is all.

        if n % 10000 == 0:
            print(f"... {n}", end="\r")
            sys.stdout.flush()

        if line[0] == 'S':
            _, name, sequence = line.strip().split('\t')

            total_bp += len(sequence)

            # first, get unitig ID
            contig_id = int(name)

            # second, make space to track the various links
            neighbors[contig_id] = []

            # third, get mean abund
            if 0:
                abund = [x for x in name_split[1:] if x.startswith("km:")]
                assert len(abund) == 1, abund
                abund = abund[0].split(":")
                assert len(abund) == 3
                abund = float(abund[2])
            abund = 0.0

            # fourth, get the min hash val for this sequence
            mh = empty_mh.copy_and_clear()
            mh.add_sequence(sequence)
            assert len(mh) == 1, (len(mh), sequence)
            min_hashval = next(iter(mh.hashes))

            # fifth, record input k-mers to a bulk signature
            in_mh.add_sequence(sequence)

            # add to database
            cursor.execute(
                "INSERT INTO sequences (id, sequence, offset, min_hashval, abund) VALUES (?, ?, ?, ?, ?)",
                (contig_id, sequence, offset, str(min_hashval), abund),
            )
        elif line[0] == 'L':
            _, fr, fro, to, too, overlap = line.strip().split()
            assert overlap == k_overlap_str # k-mer size - 1
            fr = int(fr)
            to = int(to)
            neighbors[fr].append(to)
            neighbors[to].append(fr)
        else:
            assert line[0] in ('H',)

    db.commit()

    unitigs_fp.close()

    print(f"...read {len(neighbors)} unitigs, {total_bp:.2e} bp.")
    sys.stdout.flush()

    cursor.execute("CREATE UNIQUE INDEX sequence_min_val ON sequences (min_hashval)")
    cursor.execute("CREATE UNIQUE INDEX offset_idx ON sequences (offset)")

    # sort contigs based on min_hashval!
    print("remapping cDBG IDs...")

    # remap everything into new coordinate space based on min_hashval ordering
    remapping = {}
    cursor.execute("SELECT id, offset FROM sequences ORDER BY min_hashval")
    for new_key, (old_key, offset,) in enumerate(cursor):
        remapping[old_key] = new_key

    # remap sequence IDs using offset as unique key
    cursor.execute("SELECT offset FROM sequences ORDER BY min_hashval")
    c2 = db.cursor()
    for new_key, (offset,) in enumerate(cursor):
        c2.execute("UPDATE sequences SET id=? WHERE offset=?", (new_key, offset))
        assert c2.rowcount == 1

    cursor.execute("CREATE UNIQUE INDEX sequence_idx ON sequences (id)")

    db.commit()

    print(f"DONE remapping {len(remapping)} contigs.")
    sys.stdout.flush()

    # parse link structure, map to new IDs.
    print(f"Remapping {len(neighbors)} neighbors...")
    new_neighbors = {}
    total_n = 0
    for old_key, links in neighbors.items():
        new_key = remapping[old_key]
        link_ids = set()
        for x in links:
            link_id = x
            new_link_id = remapping[link_id]
            link_ids.add(new_link_id)

        # allow isolated nodes - link_ids can be empty.
        new_neighbors[new_key] = link_ids
        total_n += len(link_ids)

    print(f"...done! {total_n} neighbor relationships.")
    sys.stdout.flush()

    # check links -- make sure that source is always in its neighbors edges.
    # (this is a check for a recurring bcalm bug that has to do with some
    # kind of threading problem)

    print("validating link structure...")
    fail = False
    for source in new_neighbors:
        for nbh in new_neighbors[source]:
            if source not in new_neighbors[nbh]:
                print(f"{source} -> {nbh}, but not {nbh} -> {source}")
                fail = True
    print("...done!")
    sys.stdout.flush()

    if fail:
        return -1

    ## save!
    print(f"saving mappings to '{args.mapping_pickle_out}'")
    with open(args.mapping_pickle_out, "wb") as fp:
        pickle.dump((ksize, new_neighbors), fp)

    # output sourmash signature for input contigs
    in_sig = sourmash.SourmashSignature(in_mh, filename=args.gfa_unitigs)
    if args.sig:
        with open(args.sig, "wt") as fp:
            sourmash.save_signatures([in_sig], fp)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
