#! /usr/bin/env python
"""
Build an index that can be used to retrieve individual reads or contigs
by cDBG node ID; produce a SQLite database for fast retrieval.

Briefly, this script creates a sqlite database with a single table,
'sequences', where a query like this:

   SELECT DISTINCT sequences.offset FROM sequences WHERE label ...

can be executed to return the offset of all sequences with the given
label. Here, 'label' is the cDBG ID to which the sequence belongs.

The script extract_reads_by_frontier_sqlite.py is a downstream script to
extract the reads with a frontier search.

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
import time
import collections
import sqlite3

import screed
from spacegraphcats.utils.logging import notify
from spacegraphcats.search import search_utils
from . import MPHF_KmerIndex, hash_sequence

# graph settings
DEFAULT_KSIZE = 31
DEFAULT_MEMORY = 1e9

def main(argv=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('reads')
    p.add_argument('savename')
    p.add_argument('-k', '--ksize', default=DEFAULT_KSIZE, type=int)
    p.add_argument('-M', '--memory', default=DEFAULT_MEMORY,
                            type=float)
    args = p.parse_args(argv)

    dbfilename = args.savename
    if os.path.exists(dbfilename):
        print('removing existing db {}'.format(dbfilename))
        os.unlink(dbfilename)

    db = sqlite3.connect(dbfilename)
    cursor = db.cursor()
    cursor.execute('CREATE TABLE sequences (offset INTEGER, cdbg_id INTEGER)');
    db.commit()

    # find the contigs filename
    contigs = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

    # ...and kmer index.
    ki_start = time.time()
    kmer_idx = MPHF_KmerIndex.from_catlas_directory(args.catlas_prefix)
    notify('loaded {} k-mers in index ({:.1f}s)',
           len(kmer_idx.table), time.time() - ki_start)

    total_bp = 0
    watermark_size = 5e5
    watermark = watermark_size

    print('walking read file: {}'.format(args.reads))
    n = 0

    cursor.execute('PRAGMA cache_size=1000000')
    cursor.execute('PRAGMA synchronous = OFF')
    cursor.execute('PRAGMA journal_mode = MEMORY')
    
    # some sqlite installs start in transactions    
    try:
        cursor.execute('BEGIN TRANSACTION')
    except sqlite3.OperationalError:
        pass

    reader = search_utils.BgzfReader(args.reads)
    total_cdbg_ids = set()
    for record, offset in search_utils.iterate_bgzf(reader):
        n += 1
        if total_bp >= watermark:
            print(f'... {watermark:5.2e} bp thru reads', end='\r',
                  file=sys.stderr)
            watermark += watermark_size
        total_bp += len(record.sequence)

        if len(record.sequence) < args.ksize:
            continue

        # identify matching cDBG IDs
        cdbg_ids = set()
        kmers = hash_sequence(record.sequence, args.ksize)
        cdbg_ids = kmer_idx.table.get_unique_values(kmers)

#        for kmer in kmers:
#            cdbg_id = kmer_idx.get_cdbg_id(kmer)
#            if cdbg_id is None:
#                continue # @CTB
#            cdbg_id = int(cdbg_id)        # convert away from numpy int
#            cdbg_ids.add(cdbg_id)
        if None in cdbg_ids: cdbg_ids.remove(None)

        for cdbg_id in cdbg_ids:
            cursor.execute('INSERT INTO sequences (offset, cdbg_id) VALUES (?, ?)', (offset, cdbg_id))

        total_cdbg_ids.update(cdbg_ids)

    db.commit()
    print(f'{watermark:5.2e} bp in reads', file=sys.stderr)
    print(f'found reads for {len(total_cdbg_ids)} cDBG IDs')
    assert max(total_cdbg_ids) + 1 == len(total_cdbg_ids)

    db.close()
    print(f'done closing {dbfilename}!')

    return 0

if __name__ == '__main__':
    sys.exit(main())
