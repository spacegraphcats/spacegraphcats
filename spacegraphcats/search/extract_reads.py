#! /usr/bin/env python
"""
Retrieve the reads for a list of cDBG nodes.  Consumes the output of
extract_nodes_by_query to get the list of nodes, and then uses the
labeled cDBG output by .cdbg.label_cdbg to find reads that overlap with
the unitigs in those nodes.
"""
import argparse
import os
import sys
import time
import gc
from collections import defaultdict
import gzip

import screed

from .search_utils import get_reads_by_cdbg

from spacegraphcats.utils.logging import log
from . import search_utils


def main(argv=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('readsfile')
    p.add_argument('labeled_reads_sqlite')
    p.add_argument('node_list_file', help='a cdbg_ids.txt.gz file')
    p.add_argument('-o', '--output')
    p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args(argv)

    dbfilename = args.labeled_reads_sqlite
    if not os.path.exists(dbfilename):
        print('sqlite file {} does not exist'.format(dbfilename))
        sys.exit(-1)


    if not args.output:
        outname = args.node_list_file + '.reads.fa.gz'
    else:
        outname = args.output

    if outname.endswith('.gz'):
        outfp = gzip.open(outname, 'wt')
    else:
        outfp = open(outname, 'wt')

    with gzip.open(args.node_list_file, 'rt') as fp:
        cdbg_shadow = set([int(x.strip()) for x in fp])

    print('extracting reads to {}.'.format(outname))

    start = time.time()

    total_bp = 0
    total_seqs = 0

    print('loading labels and querying for matching sequences..')
    reads_iter = get_reads_by_cdbg(dbfilename, args.readsfile, cdbg_shadow)
    for n, (record, offset_f) in enumerate(reads_iter):
        if n % 10000 == 0:
            print('...at n {} ({:.1f}% of file)'.format(n, offset_f * 100),
                  end='\r')

        # output!
        outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
        total_bp += len(record.sequence)
        total_seqs += 1

    print('')
    print('fetched {} reads, {} bp matching nodes.'.format(total_seqs,
                                                           total_bp))
    end = time.time()
    print('total read retrieval time (including database query): {:.2f}s'.format(end - start))

    return 0


if __name__ == '__main__':
    sys.exit(main())
