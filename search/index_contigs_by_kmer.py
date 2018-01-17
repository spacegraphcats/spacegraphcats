import sys
import sqlite3
import struct
import screed
import khmer
from .search_utils import SqlKmerIndex
import argparse
import os
from pickle import dump


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix')
    a = p.parse_args()
    
    kh = khmer.Nodetable(31, 1, 1)

    contigs_filename = os.path.join(a.catlas_prefix, 'contigs.fa.gz')
    pickle_file = os.path.join(a.catlas_prefix, 'contigs.fa.gz.kmeridx')

    d = {}
    sizes = {}

    print('reading cDBG nodes from {}'.format(contigs_filename))
    for n, record in enumerate(screed.open(contigs_filename)):
        if n % 10000 == 0 and n:
            print('...', n)

        kmers = kh.get_kmer_hashes(record.sequence)
        cdbg_id = int(record.name)

        for hashval in kmers:
            d[hashval] = cdbg_id

        sizes[cdbg_id] = len(kmers)

    with open(pickle_file, 'wb') as fp:
        dump((d, sizes), fp)


if __name__ == '__main__':
    main()
