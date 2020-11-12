#! /usr/bin/env python
"""
Produce a pickled dictionary that contains { hashval: cDBG id }, for some
specified sourmash ksize & scaled.
"""
import sys
import argparse
from pickle import dump
import sqlite3

import sourmash
import screed

from spacegraphcats.utils.logging import notify
from spacegraphcats.search import search_utils


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("contigs_db")
    parser.add_argument("picklefile")
    parser.add_argument("-k", "--ksize", type=int, default=31)
    parser.add_argument("--scaled", type=int, default=10000)
    args = parser.parse_args(argv)

    mh = sourmash.MinHash(0, args.ksize, scaled=args.scaled)
    hashval_to_contig_id = {}

    notify(f"reading contigs from {args.contigs_db}")
    contigs_db = sqlite3.connect(args.contigs_db)

    for record in search_utils.contigs_iter_sqlite(contigs_db):
        contig_id = int(record.name)

        this_mh = mh.copy_and_clear()
        this_mh.add_sequence(record.sequence, force=True)

        for hashval in this_mh.hashes:
            hashval_to_contig_id[hashval] = contig_id

    notify(
        "saving {} hashval -> cdbg_id mappings to {}",
        len(hashval_to_contig_id),
        args.picklefile,
    )
    with open(args.picklefile, "wb") as dumpfp:
        dump(hashval_to_contig_id, dumpfp)

    return 0


if __name__ == "__main__":
    main(sys.argv[1:])
