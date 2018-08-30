#! /usr/bin/env python
"""
Subtract all input sequences that have greater than some threshold overlap
with a collection of other sequences.

e.g.

    python -m search.make_donut --query nbhd1.fa nbhd2.fa \
              --subtract genome1.fa genome2.fa

would subtract genome1.fa and genome2.fa from neighborhood.fa
"""
import sys
sys.path.insert(0, '../pybbhash')
import os
from bbhash_table import BBHashTable
import khmer
import argparse
import screed


DEFAULT_THRESHOLD=0.1


def main(argv=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('--query', nargs='+', action='append')
    p.add_argument('--subtract', nargs='+', action='append')
    p.add_argument('-o', '--output-suffix')
    p.add_argument('--threshold', type=float, default=DEFAULT_THRESHOLD)
    p.add_argument('-k', '--ksize', type=int, default=31)
    p.add_argument('--directory', default='')
    args = p.parse_args(argv)

    if not args.query:
        print('error, must specify at least one query with --query')
        sys.exit(-1)

    if not args.subtract:
        print('error, must specify at least one subtract with --subtract')
        sys.exit(-1)

    args.query = [item for sublist in args.query for item in sublist]
    args.subtract = [item for sublist in args.subtract for item in sublist]

    # construct output filename as {query}.suffix
    output_suffix = args.output_suffix
    if not output_suffix:
        output_suffix = '.donut.fa'

    # load k-mers to subtract
    all_kmers = list()
    kh = khmer.Nodetable(args.ksize, 1, 1)
    
    for subtract_fn in args.subtract:
        print('loading:', subtract_fn)
        for record in screed.open(subtract_fn):
            all_kmers.extend(kh.get_kmer_hashes(record.sequence))

    # now build a minimal perfect hash function for all those k-mers
    print('building bbhash table')
    table = BBHashTable(all_kmers, fill=1)
    del all_kmers

    # next, iterate over each input and do subtract
    for queryfile in args.query:
        output = os.path.basename(queryfile) + output_suffix
        if args.directory:
            output = os.path.join(args.directory, output)
        print('subtracting from {} -> {}'.format(queryfile, output))
        outfp = open(output, 'wt')
        n = 0
        bp = 0
        n_kept = 0
        bp_kept = 0
        for n, record in enumerate(screed.open(queryfile)):
            if n % 100000 == 0:
                print('...', queryfile, n, n_kept)

            bp += len(record.sequence)

            if len(record.sequence) < args.ksize:
                continue

            kmers = kh.get_kmer_hashes(record.sequence)

            present = 0
            for k in kmers:
                if table[k]:
                    present += 1

            f = present / len(kmers)
            if f < args.threshold:            # keep?
                outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
                n_kept += 1
                bp_kept += len(record.sequence)

        print('kept {} ({:.1g} Mbp) of {} ({:.1g} Mbp)'.format(n_kept,
                                                               bp_kept / 1e6,
                                                               n, bp / 1e6))

    return 0


if __name__ == '__main__':
    sys.exit(main())
