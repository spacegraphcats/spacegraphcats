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
from . import bgzf

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
    cursor.execute('CREATE TABLE sequences (offset INTEGER, label INTEGER)');
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
    tags_to_label = collections.defaultdict(int)
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
            tags_to_label[t] = cdbg_id

    ###

    total_bp = 0
    watermark_size = 1e7
    watermark = watermark_size

    def read_bgzf(filename):
        from screed.openscreed import fastq_iter, fasta_iter
        reader = bgzf.BgzfReader(filename, 'rt')
        ch = reader.read(1)
        if ch == '>':
            iter_fn = fasta_iter
        elif ch == '@':
            iter_fn = fastq_iter
        else:
            raise Exception('unknown start chr {}'.format(ch))

        reader.seek(0)

        last_pos = reader.tell()
        for record in iter_fn(reader):
            yield record, last_pos
            last_pos = reader.tell()

    print('walking read file: {}'.format(args.reads))
    n = 0

    cursor.execute('PRAGMA cache_size=1000000')
    cursor.execute('PRAGMA synchronous = OFF')
    cursor.execute('PRAGMA journal_mode = MEMORY')
    cursor.execute('BEGIN TRANSACTION')

    for record, offset in read_bgzf(args.reads):
        n += 1
        if total_bp >= watermark:
            print('... {:5.2e} bp thru reads'.format(int(watermark)),
                  file=sys.stderr)
            watermark += watermark_size
        total_bp += len(record.sequence)

        tags = ng.get_tags_for_sequence(record.sequence)
        labels = set([ tags_to_label[t] for t in tags ])

        for lb in labels:
            cursor.execute('INSERT INTO sequences (offset, label) VALUES (?, ?)', (offset, lb))

    db.commit()
            
    db.close()

if __name__ == '__main__':
    main()
