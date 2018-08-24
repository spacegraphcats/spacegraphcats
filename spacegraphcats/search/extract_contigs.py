#! /usr/bin/env python
"""
Retrieve the unitig sequences for a given list of cDBG nodes.  Consumes the
output of extract_nodes_by_query to get the list of nodes.
"""
import argparse
import os
import sys
import gzip

import screed

from spacegraphcats.utils.logging import log
from . import search_utils


def main(argv=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('node_list_file', help='a cdbg_ids.txt.gz file')
    p.add_argument('-o', '--output')
    p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args(argv)

    contigs = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

    if not args.output:
        outname = args.node_list_file + '.contigs.fa.gz'
    else:
        outname = args.output

    if outname.endswith('.gz'):
        outfp = gzip.open(outname, 'wt')
    else:
        outfp = open(outname, 'wt')

    with gzip.open(args.node_list_file, 'rt') as fp:
        cdbg_shadow = set([ int(x.strip()) for x in fp ])

    total_bp = 0
    total_seqs = 0

    print('extracting contigs to {}.'.format(outname))
    for n, record in enumerate(screed.open(contigs)):
        if n % 10000 == 0:
            offset_f = total_seqs / len(cdbg_shadow)
            print('...at n {} ({:.1f}% of shadow)'.format(total_seqs,
                                                          offset_f * 100),
                  end='\r')

        contig_id = int(record.name)
        if contig_id not in cdbg_shadow:
            continue

        outfp.write('>{}\n{}\n'.format(record.name, record.sequence))

        total_bp += len(record.sequence)
        total_seqs += 1

    print('')
    print('fetched {} contigs, {} bp matching node list.'.format(total_seqs, total_bp))

    return 0


if __name__ == '__main__':
    sys.exit(main())
