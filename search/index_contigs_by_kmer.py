import sys
import sqlite3
import struct
import screed
import khmer
from .search_utils import SqlKmerIndex
import argparse
import os


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix')
    a = p.parse_args()
    
    kh = khmer.Nodetable(31, 1, 1)

    contigs_filename = os.path.join(a.catlas_prefix, 'contigs.fa.gz')

    db = SqlKmerIndex(contigs_filename + '.sqlindex')
    db.create()

    print('reading cDBG nodes from {}'.format(contigs_filename))
    for n, record in enumerate(screed.open(contigs_filename)):
        if n % 10000 == 0 and n:
            print('...', n)

        kmers = kh.get_kmer_hashes(record.sequence)
        cdbg_id = int(record.name)

        for hashval in kmers:
            db.insert_kmer(hashval, cdbg_id)

        db.set_node_size(cdbg_id, len(kmers))

    print('committing...')
    db.commit()
    print('...building index...')
    db.build_index()
    print('...committing')
    db.commit()


if __name__ == '__main__':
    main()
