#! /usr/bin/env python
import sys
import os
import argparse
import screed
import khmer
import collections
from pickle import dump

sys.path.insert(0, '/Users/t/dev/khmer')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('savename')
    p.add_argument('-k', '--ksize', default=31, type=int)
    args = p.parse_args()

    ng = khmer.Nodegraph(args.ksize, 2e9, 4)

    basename = os.path.basename(args.catlas_prefix)
    contigfile = os.path.join(args.catlas_prefix, "contigs.fa.gz")

    total_bp = 0
    watermark_size = 1e6
    watermark = watermark_size

    tags_to_labels = collections.defaultdict(set)

    n = 0
    for record in screed.open(contigfile):
        n += 1
        if total_bp >= watermark:
            print('... {:5.2e} bp thru contigs'.format(int(watermark)),
                  file=sys.stderr)
            watermark += watermark_size
        total_bp += len(record.sequence)

        cdbg_id = int(record.name)
        ng.consume_and_tag(record.sequence)
        for t in ng.get_tags_for_sequence(record.sequence):
            tags_to_labels[t].add(cdbg_id)
        
        if 0:
            print('----')
            print('bar', record.name, len(record.sequence))
            print('foo',
                  len(lh.get_tags_for_sequence(record.sequence)),
                  len(lh.get_labels_for_sequence(record.sequence)),
                  lh.get_labels_for_sequence(record.sequence),
                  cdbg_id)
            print('foo2', list(lh.get_tags_for_sequence(record.sequence)))

    ng.save(args.savename)
    ng.save_tagset(args.savename + '.tagset')
    dump(tags_to_labels, open(args.savename + '.labelsp', 'wb'))


if __name__ == '__main__':
    main()
