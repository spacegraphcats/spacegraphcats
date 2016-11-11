#! /usr/bin/env python
"""
Take a list of murmurhashes output by gimme-dbg-nodes and sweep through
a file of reads, extracting all of the reads that contain k-mers with
those hashes into a separate output file.
"""
import khmer, khmer.utils
from sourmash_lib._minhash import hash_murmur
import argparse
import os
import screed


K=31


def kmers(seq, ksize):
    seq = seq.upper()
    for start in range(len(seq) - ksize + 1):
        yield seq[start:start + ksize]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('readfile')
    parser.add_argument('hashfile')
    parser.add_argument('-o', '--output', metavar="filename",
                        type=argparse.FileType('wb'),
                        default=None, dest='output')
    args = parser.parse_args()

    with open(args.hashfile, 'r') as fp:
        hashes = set()
        for line in fp:
            hashes.add(int(line))

    if args.output:
        output_fp = args.output
    else:
        output_filename = os.path.basename(args.readfile) + '.sweep'
        output_fp = open(output_filename, 'w')

    m = 0
    for n, record in enumerate(screed.open(args.readfile)):
        keep = False
        for k in kmers(record.sequence, K):
            if hash_murmur(k) in hashes:
                keep = True
                break

        if keep:
            m += 1
            khmer.utils.write_record(record, output_fp)

    print('kept %d of %d reads, for %d hashes' % (m, n, len(hashes)))


if __name__ == '__main__':
    main()
