#! /usr/bin/env python
"""
Dump the data in the contigs sqlite db to FASTA.
"""
import sys
import argparse
import sqlite3
from spacegraphcats.search import search_utils


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("sqlite_db")
    args = parser.parse_args(argv)

    db = sqlite3.connect(args.sqlite_db)
    for record in search_utils.contigs_iter_sqlite(db):
        print(f">{record.name}\n{record.sequence}")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
