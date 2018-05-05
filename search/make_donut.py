#! /usr/bin/env python
"""
Subtract all input sequences that have greater than some threshold overlap
with a collection of other sequences.

e.g.

    python -m search.make_donut neighborhood.fa genome1.fa genome2.fa

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
    p.add_argument('sequences')
    p.add_argument('subtract', nargs='+')
    p.add_argument('-o', '--output')
    p.add_argument('--threshold', type=float, default=DEFAULT_THRESHOLD)
    p.add_argument('-k', '--ksize', type=int, default=31)
    args = p.parse_args(argv)

    # construct output filename as sequences.donut.fa
    output = args.output
    if not output:
        output = os.path.basename(args.sequences) + '.donut.fa'

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

    # next, iterate over input and do subtract
    print('subtracting from {} -> {}'.format(args.sequences, output))
    outfp = open(output, 'wt')
    n = 0
    bp = 0
    n_kept = 0
    bp_kept = 0
    for n, record in enumerate(screed.open(args.sequences)):
        if n % 100000 == 0:
            print('...', n, n_kept)

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
