#! /usr/bin/env python3

from khmer import MinHash
import screed
import os.path
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqfile')
    parser.add_argument('-o', '--out',
                        metavar='outfile', type=argparse.FileType('wt'))
    args = parser.parse_args()

    mh = MinHash(1000, 31)

    print('loading sequences from', args.seqfile)
    for record in screed.open(args.seqfile):
        mh.add_sequence(record.sequence)

    if not args.out:
        outfilename = os.path.basename(args.seqfile) + '.dump.txt'
        outfile = open(outfilename, 'wt')
    else:
        outfile = args.out

    print('dumping sequences to', outfile.name)
    outfile.write(' '.join(map(str, mh.get_mins())))
    outfile.close()
    print('done!')

        
if __name__ == '__main__':
    main()
