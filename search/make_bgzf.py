#! /usr/bin/env python
"""
Convert an input file into a BGZF file that supports random indexing by offset.
"""
import screed
from search.bgzf import bgzf
import os.path
import argparse


def main():
    p = argparse.ArgumentParser()
    p.add_argument('input_file')
    args = p.parse_args()

    output_filename = os.path.basename(args.input_file) + '.bgz'
    outfp = bgzf.open(output_filename, 'wb')

    print('turning {} into a block-gzipped (BGZF) file'.format(args.input_file))
    print('output file will be {}'.format(output_filename))

    for n, record in enumerate(screed.open(args.input_file)):
        offset = outfp.tell()
        outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
        if n % 10000 == 0:
            print('offset for {} is {}'.format(n, offset))

    outfp.close()

    print('done!')


if __name__ == '__main__':
    main()
