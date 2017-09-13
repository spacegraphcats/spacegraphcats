#! /usr/bin/env python
"""
Build an index that can be used to retrieve individual reads or contigs
by cDBG node ID; produce a SQLite database for fast retrieval.

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
import sqlite3
import shutil

sys.path.insert(0, '/Users/t/dev/khmer')

# graph settings
DEFAULT_KSIZE = 31
DEFAULT_MEMORY = 1e9

def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('reads')
    p.add_argument('savename')
    p.add_argument('-k', '--ksize', default=DEFAULT_KSIZE, type=int)
    p.add_argument('-M', '--memory', default=DEFAULT_MEMORY,
                            type=float)
    args = p.parse_args()

    dbfilename = args.savename + '.sqlite'
    if os.path.exists(dbfilename):
        print('removing existing db {}'.format(dbfilename))
        os.unlink(dbfilename)

    db = sqlite3.connect(dbfilename)
    cursor = db.cursor()
    cursor.execute('CREATE TABLE sequences (id INTEGER PRIMARY KEY, name TEXT, sequence TEXT, quality TEXT)');
    cursor.execute('CREATE TABLE tags_to_sequences (tag INTEGER, seq_id INTEGER, FOREIGN KEY(seq_id) REFERENCES sequences(id))');
    cursor.execute('CREATE TABLE tags_and_labels(tag INTEGER, label INTEGER)');
    db.commit()

    # @CTB support different sizes.
    graph_tablesize = int(args.memory * 8.0 / 4.0)
    ng = khmer.Nodegraph(args.ksize, graph_tablesize, 4)

    basename = os.path.basename(args.catlas_prefix)
    contigfile = os.path.join(args.catlas_prefix, "contigs.fa.gz")

    total_bp = 0
    watermark_size = 1e6
    watermark = watermark_size

    print('walking catlas cDBG contigs: {}'.format(contigfile))
    n = 0
    for contig in screed.open(contigfile):
        n += 1
        if total_bp >= watermark:
            print('... {:5.2e} bp thru contigs'.format(int(watermark)),
                  file=sys.stderr)
            watermark += watermark_size
        total_bp += len(contig.sequence)

        cdbg_id = int(contig.name)
        ng.consume_and_tag(contig.sequence)

        tags = ng.get_tags_for_sequence(contig.sequence)

        for t in tags:
            cursor.execute('INSERT INTO tags_and_labels (tag, label) VALUES (?, ?)', (t, cdbg_id))

    db.commit()

    ###

    total_bp = 0
    watermark_size = 1e7
    watermark = watermark_size

    print('walking read file: {}'.format(args.reads))
    n = 0
    for record in screed.open(args.reads):
        n += 1
        if total_bp >= watermark:
            print('... {:5.2e} bp thru reads'.format(int(watermark)),
                  file=sys.stderr)
            watermark += watermark_size
        total_bp += len(record.sequence)

        if not hasattr(record, 'quality'):
            record.quality = ''

        cursor.execute('INSERT INTO sequences (id, name, sequence, quality) VALUES (?, ?, ?, ?)', (n, record.name, record.sequence, record.quality))

        tags = ng.get_tags_for_sequence(record.sequence)
        for t in tags:
            cursor.execute('INSERT INTO tags_to_sequences (tag, seq_id) VALUES (?, ?)', (t, n))

    db.commit()
            
    db.close()

if __name__ == '__main__':
    main()
