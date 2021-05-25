"""
Use Minimal Perfect Hashing (see BBHash) to construct a fast lookup
table connecting k-mers in the cDBG to cDBG node IDs.

Input: a SQLite database containing the contigs in question.

Uses the bbhash package to construct a table mapping k-mer -> contigID,
and saves that as as well as a pickled file full of contig sizes.

Note: relies on the fact that for a cDBG constructed at a particular k,
no k-mer will appear in more than one cDBG node and every k-mer will
be in at least one cDBG node, i.e. k <-> cdbg_id is bijective.
"""
import sys
import argparse
import pickle
import sqlite3

import screed
from bbhash_table import BBHashTable

from .hashing import hash_sequence, MPHF_KmerIndex

UNSET_VALUE = 2**32 - 1             # BBHashTable "unset" value


def build_mphf(ksize, records_iter_fn):
    # build a list of all k-mers in the cDBG
    all_kmers = set()
    sum_kmers = 0
    multicounts = set()

    records_iter = records_iter_fn()
    for n, record in enumerate(records_iter):
        if n % 50000 == 0 and n:
            print("... contig", n, end="\r")

        kmers = hash_sequence(record.sequence, ksize)
        sum_kmers += len(kmers)
        these_kmers = set(kmers)

        # have we seen any of these kmers before? if so, remove.
        seen_before = all_kmers.intersection(these_kmers)
        if seen_before:
            multicounts.update(seen_before)

        all_kmers.update(these_kmers)

    n_contigs = n + 1
    print(f"loaded {n_contigs} contigs.\n")

    if multicounts:
        print('NOTE: likely hash collisions (or duplicate k-mers?) in input cDBG')
        print(f'NOTE: {len(multicounts)} k-mer hash values are present more than once.')
        print('NOTE: these k-mers are being removed from consideration.')

        all_kmers -= multicounts
    else:
        print('NOTE: no multicount hashvals detected.')

    # build MPHF (this is the CPU intensive bit)
    print(f"building MPHF for {len(all_kmers)} k-mers in {n_contigs} nodes.")
    table = BBHashTable()
    table.initialize(list(all_kmers))

    del all_kmers

    # build tables linking:
    # * mphf hash to k-mer hash (for checking exactness)
    # * mphf hash to cDBG ID
    # * cDBG ID to node size (in k-mers)

    print("second pass.")
    records_iter = records_iter_fn()
    sizes = {}
    max_cdbg_id = 0
    for n, record in enumerate(records_iter):
        if n % 50000 == 0 and n:
            print("... contig {} of {}".format(n, n_contigs), end="\r")

        # node ID is record name, must go from 0 to total-1
        cdbg_id = int(record.name)

        # get 64-bit numbers for each k-mer
        kmers = hash_sequence(record.sequence, ksize)

        # for each k-mer, find its MPHF hashval, & link to info.
        for kmer in set(kmers) - multicounts:
            table[kmer] = cdbg_id

        sizes[cdbg_id] = len(kmers)

        # update max cdbg_id:
        if cdbg_id > max_cdbg_id:
            max_cdbg_id = cdbg_id

    print(f"loaded {n_contigs} contigs in pass2.\n")

    # a few asserts - these are a bit redundant with each other, but
    # are here to aid in debugging.

    # number of contigs should not be super large.
    assert n <= UNSET_VALUE

    # no values in the bbhash table should be unset.
    max_table_value = max(table.mphf_to_value)
    assert max_table_value != UNSET_VALUE

    # cdbg_id should go from 0 to n_contigs - 1.
    assert n == max_cdbg_id

    # and max value in the bbhash table should be the max cdbg_id
    assert n == max_table_value, (n, max_table_value)

    return table, sizes


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument("catlas_prefix")
    p.add_argument("--contigs-db", required=True)
    p.add_argument("-k", "--ksize", default=31, type=int)
    args = p.parse_args(argv)

    sqlite_db = sqlite3.connect(args.contigs_db)

    def create_records_iter():
        from spacegraphcats.search import search_utils

        print(f"reading cDBG nodes from sqlite DB {args.contigs_db}")
        return search_utils.contigs_iter_sqlite(sqlite_db)

    table, sizes = build_mphf(args.ksize, create_records_iter)
    print(f"done! saving to {args.catlas_prefix}")

    kmer_idx = MPHF_KmerIndex(args.ksize, table, sizes)
    kmer_idx.save(args.catlas_prefix)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
