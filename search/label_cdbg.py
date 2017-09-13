#! /usr/bin/env python
"""
Build an index that can be used to retrieve individual reads or contigs
by cDBG node ID.

Specifically,
* walk through the contigs assembled from the cDBG;
* build a DBG cover using khmer tags, such that every k-mer in the DBG
  is within distance d=40 of a tag;
* label each tag with the cDBG node ID from the contig;
* save for later use.
"""
import sys
import os
import argparse
import screed
import khmer
import collections
from pickle import dump

sys.path.insert(0, '/Users/t/dev/khmer')

# graph settings
DEFAULT_KSIZE = 31
DEFAULT_MEMORY = 1e9

def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('savename')
    p.add_argument('-k', '--ksize', default=DEFAULT_KSIZE, type=int)
    p.add_argument('-M', '--memory', default=DEFAULT_MEMORY,
                            type=float)
    args = p.parse_args()

    graph_tablesize = int(args.memory * 8.0 / 4.0)
    ng = khmer.Nodegraph(args.ksize, graph_tablesize, 4)

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

        tags = ng.get_tags_for_sequence(record.sequence)
        for t in tags:
            tags_to_labels[t].add(cdbg_id)

    ng.save_tagset(args.savename + '.tagset')
    dump(tags_to_labels, open(args.savename + '.labelsp', 'wb'))

    print('saved {} tags, {} labels for {} bp'.format(ng.n_tags(), n, total_bp))


if __name__ == '__main__':
    main()
