#! /usr/bin/env python
"""
Build an index that can be used to retrieve individual reads or contigs
by cDBG node ID; produce a SQLite database for fast retrieval.

Briefly, this script creates a sqlite database with a single table,
'sequences', where a query like this:

   SELECT DISTINCT sequences.offset FROM sequences WHERE label ...

can be executed to return the offset of all sequences with the given
label. Here, 'label' is the cDBG ID to which the sequence belongs.

The script query_by_sequence is a downstream script to extract the
reads.

Specifically,
* walk through the reads and find their offsets in the BGZF file
* use the k-mer index created by `index_cdbg_by_kmer` to assign them to
  contig(s)
* record offset of reads <-> cDBG ID
"""
import sys
import os
import argparse
import time
import collections
import sqlite3

import screed
from spacegraphcats.utils.logging import notify
from spacegraphcats.utils.bgzf.bgzf import BgzfReader
from spacegraphcats.search import search_utils
from . import MPHF_KmerIndex, hash_sequence

# graph settings
DEFAULT_KSIZE = 31


def main(argv=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("reads")
    p.add_argument("savename")
    p.add_argument("-k", "--ksize", default=DEFAULT_KSIZE, type=int)
    p.add_argument("-P", "--expect-paired", action="store_true")
    p.add_argument("-N", "--ignore-paired", action="store_true")
    args = p.parse_args(argv)

    if args.expect_paired and args.ignore_paired:
        print("cannot set both -P/--expect-paired and -N/--ignore-paired")
        return -1

    if args.expect_paired:
        print(
            "We will REQUIRE that some of the reads are in pairs (-P/--expect-paired)"
        )
    else:
        print("We will NOT require that some of theads be in pairs (default).")

    if args.ignore_paired:
        print("Ignoring paired reads (-N/--ignore-paired)")
    else:
        print("Paired reads will be indexed together (default).")

    dbfilename = args.savename
    if os.path.exists(dbfilename):
        print(f"removing existing db '{dbfilename}'")
        os.unlink(dbfilename)

    db = sqlite3.connect(dbfilename)
    cursor = db.cursor()
    cursor.execute("CREATE TABLE sequences (offset INTEGER, cdbg_id INTEGER)")
    db.commit()

    # load k-mer MPHF index: { kmers -> cDBG IDs }
    ki_start = time.time()
    kmer_idx = MPHF_KmerIndex.from_catlas_directory(args.catlas_prefix)
    notify("loaded {} k-mers in index ({:.1f}s)", len(kmer_idx), time.time() - ki_start)

    total_bp = 0
    watermark_size = 5e5
    watermark = watermark_size

    print(f"walking read file: '{args.reads}'")
    n = 0

    cursor.execute("PRAGMA cache_size=1000000")
    cursor.execute("PRAGMA synchronous = OFF")
    cursor.execute("PRAGMA journal_mode = MEMORY")

    # some sqlite installs start in transactions
    try:
        cursor.execute("BEGIN TRANSACTION")
    except sqlite3.OperationalError:
        pass

    reader = BgzfReader(args.reads)
    total_cdbg_ids = set()
    last_record = None
    last_offset = None
    n_paired_reads = 0
    for record, offset in search_utils.iterate_bgzf(reader):
        n += 1
        if total_bp >= watermark:
            print(f"... {watermark:5.2e} bp thru reads", end="\r", file=sys.stderr)
            watermark += watermark_size

            if args.expect_paired:
                if not n_paired_reads:
                    print(
                        "ERROR: no paired reads!? but -P/--expect-paired set. Quitting."
                    )
                    return -1

        total_bp += len(record.sequence)

        this_name = record.name
        is_paired = False
        if last_record:
            last_name = last_record.name

            if last_name.endswith("/1") and this_name.endswith("/2"):
                if last_name[:-2] == this_name[:-2] and len(last_name) > 2:
                    is_paired = True
            elif last_name == this_name and last_name:  # SRR stuff, bleah
                is_paired = True

        if is_paired:
            offsets = [last_offset, offset]
            # leave cdbg_ids alone... will be cleared next.
            n_paired_reads += 1
        else:
            offsets = [offset]
            cdbg_ids = set()

        if len(record.sequence) < args.ksize:
            continue

        # identify matching cDBG IDs
        kmers = hash_sequence(record.sequence, args.ksize)

        # CTB note: if reads have 'N' in them, they will not appear
        # in cDBG, so do not require_exist=True here.
        new_cdbg_ids = kmer_idx.count_cdbg_matches(kmers)
        cdbg_ids.update(new_cdbg_ids)

        for insert_offset in offsets:
            # CTB: note this may insert duplicates into the table,
            # which are resolved in the extraction by SELECT DISTINCT.
            for cdbg_id in cdbg_ids:
                cursor.execute(
                    "INSERT INTO sequences (offset, cdbg_id) VALUES (?, ?)",
                    (insert_offset, cdbg_id),
                )

        total_cdbg_ids.update(cdbg_ids)

        if not args.ignore_paired:
            last_record = record
            last_offset = offset

    db.commit()
    notify(f"{total_bp:5.2e} bp in {n} reads")
    notify(f"{n_paired_reads} paired of {n} reads")
    notify(f"found reads for {len(total_cdbg_ids)} cDBG IDs")
    assert max(total_cdbg_ids) + 1 == len(total_cdbg_ids)

    cursor.execute("CREATE INDEX offset_idx ON sequences (offset)")
    cursor.execute("CREATE INDEX cdbg_id_idx ON sequences (cdbg_id)")
    db.close()
    print(f"done; closing {dbfilename}!")

    if args.expect_paired:
        if not n_paired_reads:
            print("ERROR: no paired reads!? but -P/--expect-paired set. Failing.")
            return -1

    return 0


if __name__ == "__main__":
    sys.exit(main())
