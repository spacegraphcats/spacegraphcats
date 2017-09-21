#! /usr/bin/env python
"""
Convert an input file into a BGZF file that supports random indexing by offset.
"""
import sys
import screed
from search.bgzf import bgzf
import os.path
import argparse


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('input_file')
    args = p.parse_args()

    outfp = bgzf.open(os.path.basename(args.input_file) + '.bgz', 'wb')

    for n, record in enumerate(screed.open(sys.argv[1])):
        offset = outfp.tell()
        outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
        if n % 10000 == 0:
            print('offset for {} is {}'.format(n, offset))


if __name__ == '__main__':
    main()
