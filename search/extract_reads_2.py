#! /usr/bin/env python
import sys
import os
import argparse
import screed
import khmer, khmer.utils
import collections
from pickle import load, dump


def main():
    p = argparse.ArgumentParser()
    p.add_argument('labelsname')
    p.add_argument('reads_inp')
    p.add_argument('reads_outp')
    p.add_argument('-k', '--ksize', default=31, type=int)
    args = p.parse_args()

    print('loading...')
    # actually don't need whole graph, just tags @CTB
    ng = khmer.load_nodegraph(args.labelsname)
    ng.load_tagset(args.labelsname + '.tagset')

    tags_to_labels = load(open(args.labelsname + '.labelsp', 'rb'))

    total_bp = 0
    watermark_size = 1e6
    watermark = watermark_size
    no_tags = 0
    no_tags_bp = 0
    total_seqs = 0

    outfp = open(args.reads_outp, 'wt')

    n = 0
    for record in screed.open(args.reads_inp):
        n += 1
        if total_bp >= watermark:
            print('... {:5.2e} bp thru reads'.format(int(watermark)),
                  file=sys.stderr)
            watermark += watermark_size

        total_bp += len(record.sequence)
        total_seqs += 1

        tags = ng.get_tags_for_sequence(record.sequence)

        labels = set()
        for t in tags:
            labels.update(tags_to_labels[t])

            
        if not len(labels):
            no_tags += 1
            no_tags_bp += len(record.sequence)

        record.name = ",".join([str(x) for x in labels])
        
        khmer.utils.write_record(record, outfp)

    print(no_tags, total_seqs)
    print('{:5.2e} {:5.2e}'.format(no_tags_bp, total_bp))


if __name__ == '__main__':
    main()
