#! /usr/bin/env python
"""
Sort the bcalm unitigs.fa output (a cDBG) into deterministic order.

Outputs a mapping rather than a new file!

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
    parser.add_argument("bcalm_unitigs")
    parser.add_argument("sqlite_db_out")
    parser.add_argument("mapping_pickle_out")
    parser.add_argument("-k", "--ksize", type=int, default=31)
    args = parser.parse_args(argv)

    unitigs = args.bcalm_unitigs
    ksize = args.ksize

    empty_mh = sourmash.MinHash(n=1, ksize=ksize)

    ###
    db = sqlite3.connect(args.sqlite_db_out)
    cursor = db.cursor()
    cursor.execute("PRAGMA synchronous='OFF'")
    cursor.execute("PRAGMA locking_mode=EXCLUSIVE")
    cursor.execute('CREATE TABLE sequences (id INTEGER, sequence TEXT, offset INTEGER PRIMARY KEY, min_hashval TEXT, abund FLOAT)')

    # read & find minimum hash:
    total_bp = 0
    neighbors = collections.defaultdict(set)

    # record input k-mers in a minhash
    in_mh = sourmash.MinHash(0, ksize, scaled=1000)

    print(f"loading BCALM unitigs from {unitigs}...")

    unitigs_fp = open(unitigs, "rt")
    for n, (record, offset) in enumerate(my_fasta_iter(unitigs_fp)):
        if n % 1000 == 0:
            print(f'... {n}', end='\r')
            sys.stdout.flush()
        total_bp += len(record.sequence)

        # first, get unitig ID
        name = record.name
        name_split = name.split()
        contig_id = int(name_split[0])

        # second, track the various links
        links = [x for x in name_split[1:] if x.startswith("L:")]
        link_ids = [x.split(":")[2] for x in links]
        link_ids = [int(x) for x in link_ids if int(x) != contig_id]

        neighbors[contig_id].update(link_ids)

        # third, get mean abund
        abund = [x for x in name_split[1:] if x.startswith("km:")]
        assert len(abund) == 1, abund
        abund = abund[0].split(":")
        assert len(abund) == 3
        abund = float(abund[2])

        # fourth, get the min hash val for this sequence
        mh = empty_mh.copy_and_clear()
        mh.add_sequence(record.sequence)
        assert len(mh) == 1, (len(mh), record.sequence)
        min_hashval = next(iter(mh.hashes))

        # fifth, record input k-mers to a bulk signature
        in_mh.add_sequence(record.sequence)

        # add to database
        cursor.execute('INSERT INTO sequences (id, sequence, offset, min_hashval, abund) VALUES (?, ?, ?, ?, ?)', (contig_id, record.sequence, offset, str(min_hashval), abund))

    db.commit()

    print(f"...read {len(neighbors)} unitigs, {total_bp:.2e} bp.")
    sys.stdout.flush()

    # check links -- make sure that source is always in its neighbors edges.
    # (this is a check for a recurring bcalm bug that has to do with some
    # kind of threading problem)
    print("validating link structure...")
    fail = False
    for source in neighbors:
        for nbhd in neighbors[source]:
            if source not in neighbors[nbhd]:
                print(f"{source} -> {nbhd}, but not {nbhd} -> {source}")
                fail = True
    print("...done!")
    sys.stdout.flush()

    cursor.execute('CREATE UNIQUE INDEX sequence_min_val ON sequences (min_hashval)')

    if fail:
        return -1

    # sort contigs based on min_hashval!
    print("remapping cDBG IDs...")

    # remap everything into new coordinate space
    remapping = {}
    cursor.execute('SELECT id, offset FROM sequences ORDER BY min_hashval')
    for new_key, (old_key, offset,) in enumerate(cursor):
        remapping[old_key] = new_key

    # remap sequence IDs
    cursor.execute('SELECT offset FROM sequences ORDER BY min_hashval')
    c2 = db.cursor()
    for new_key, (offset,) in enumerate(cursor):
        c2.execute('UPDATE sequences SET id=? WHERE offset=?',
                   (new_key, offset))
        assert c2.rowcount == 1

    cursor.execute('CREATE UNIQUE INDEX sequence_idx ON sequences (id)')

    db.commit()

    print(f'DONE remapping {len(remapping)} contigs.')
    sys.stdout.flush()

    # remap other things
    print(f'Remapping {len(neighbors)} neighbors...')
    new_neighbors = collections.defaultdict(set)
    for old_key, vv in neighbors.items():
        new_vv = [remapping[v] for v in vv]
        new_neighbors[remapping[old_key]] = set(new_vv)
    neighbors = new_neighbors

    print("...done!")
    sys.stdout.flush()

    ## save!
    print(f"saving mappings to '{args.mapping_pickle_out}'")
    with open(args.mapping_pickle_out, "wb") as fp:
        pickle.dump((ksize, neighbors), fp)

    # output sourmash signature for input contigs
    in_sig = sourmash.SourmashSignature(in_mh, filename=args.bcalm_unitigs)
    sourmash.save_signatures([in_sig], open(args.bcalm_unitigs + ".sig", "wt"))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
