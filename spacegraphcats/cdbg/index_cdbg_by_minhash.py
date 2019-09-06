#! /usr/bin/env python
"""
Produce a pickled dictionary that contains hashval: cDBG id, for some
specified sourmash ksize & scaled.
"""
import screed
import sys
import collections
import argparse
import sourmash
from pickle import dump


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('contigs')
    parser.add_argument('picklefile')
    parser.add_argument('-k', '--ksize', type=int, default=31)
    parser.add_argument('--scaled', type=int, default=10000)
    args = parser.parse_args(argv)

    mh = sourmash.MinHash(0, args.ksize, scaled=args.scaled)
    hashval_to_contig_id = {}

    for record in screed.open(args.contigs):
        contig_id = int(record.name)

        this_mh = mh.copy_and_clear()
        this_mh.add_sequence(record.sequence, force=True)
        mins = this_mh.get_mins()

        for hashval in mins:
            hashval_to_contig_id[hashval] = contig_id

    with open(args.picklefile, 'wb') as dumpfp:
        dump(hashval_to_contig_id, dumpfp)


if __name__ == '__main__':
    main(sys.argv[1:])
