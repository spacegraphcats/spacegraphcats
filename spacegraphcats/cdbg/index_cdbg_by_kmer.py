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
import os
import argparse
import pickle
import sqlite3

import screed
import khmer
from bbhash_table import BBHashTable


hashing_fn = None
hashing_ksize = None


def hash_sequence(seq, ksize):
    global hashing_fn, hashing_ksize
    if hashing_fn is None or hashing_ksize != ksize:
        kh = khmer.Nodetable(ksize, 1, 1)
        hashing_fn, hashing_ksize = kh.get_kmer_hashes, ksize

    return hashing_fn(seq)


class MPHF_KmerIndex(object):
    """
    Support kmer -> cDBG node id queries.
    """

    def __init__(self, bbhash_table, sizes):
        self.table = bbhash_table
        self.sizes = sizes

    def __len__(self):
        return len(self.table)

    def get_cdbg_id(self, kmer_hash):
        return self.table[kmer_hash]

    def get_cdbg_size(self, cdbg_id):
        return self.sizes[cdbg_id]

    def build_catlas_node_sizes(self, catlas):
        """Count the number of kmers in the shadow of each catlas node."""
        node_kmer_sizes = {}

        # (this iteration is always in order of levels, note)
        for node_id in catlas:
            level = catlas.levels[node_id]

            # for level 1 nodes, get the total k-mers in cDBG nodes
            # under this dominator.
            if level == 1:
                total_kmers = 0
                for cdbg_node in catlas.layer1_to_cdbg.get(node_id):
                    total_kmers += self.get_cdbg_size(cdbg_node)

                node_kmer_sizes[node_id] = total_kmers
            # for higher level nodes, it's the sum of sizes of the children
            else:
                sub_size = 0
                for child_id in catlas.children[node_id]:
                    sub_size += node_kmer_sizes[child_id]
                node_kmer_sizes[node_id] = sub_size

        return node_kmer_sizes

    def count_cdbg_matches(self, query_kmers, verbose=True, require_exist=False):
        "Return a dictionary containing cdbg_id -> # of matches in query_kmers"
        return self.table.get_unique_values(query_kmers, require_exist=require_exist)

    def count_catlas_matches(self, cdbg_matches, catlas):
        """ """
        dom_matches = {}
        total_kmers_in_query_nodes = 0.0

        # for every node in the catlas,
        for node_id in catlas:
            level = catlas.levels[node_id]

            # if bottom dominator, get the counts for k-mers matched
            # to all the cDBG nodes under it.
            if level == 1:
                query_kmers = 0
                all_kmers_in_node = 0

                # sum both matching k-mers and all k-mers in cDBG nodes
                for cdbg_node in catlas.layer1_to_cdbg.get(node_id):
                    query_kmers += cdbg_matches.get(cdbg_node, 0)
                    all_kmers_in_node += self.get_cdbg_size(cdbg_node)

                # if anything found, save -> dom_matches
                if query_kmers:
                    dom_matches[node_id] = query_kmers
                    total_kmers_in_query_nodes += all_kmers_in_node

                    assert all_kmers_in_node >= query_kmers

            else:
                # propogate up the catlas.
                sub_size = 0
                for child_id in catlas.children[node_id]:
                    sub_size += dom_matches.get(child_id, 0)

                if sub_size:
                    dom_matches[node_id] = sub_size

        return dom_matches

    @classmethod
    def from_catlas_directory(cls, catlas_prefix):
        "Load kmer index created by search.contigs_by_kmer."
        mphf_filename = os.path.join(catlas_prefix, "contigs.mphf")
        array_filename = os.path.join(catlas_prefix, "contigs.indices")
        sizes_filename = os.path.join(catlas_prefix, "contigs.sizes")

        table = BBHashTable.load(mphf_filename, array_filename)
        with open(sizes_filename, "rb") as fp:
            sizes = pickle.load(fp)

        return cls(table, sizes)


def build_mphf(ksize, records_iter_fn):
    # build a list of all k-mers in the cDBG
    all_kmers = list()

    records_iter = records_iter_fn()
    for n, record in enumerate(records_iter):
        if n % 50000 == 0 and n:
            print("... contig", n, end="\r")

        kmers = hash_sequence(record.sequence, ksize)
        all_kmers.extend(list(kmers))

    n_contigs = n + 1
    print(f"loaded {n_contigs} contigs.\n")

    # build MPHF (this is the CPU intensive bit)
    print(f"building MPHF for {len(all_kmers)} k-mers in {n_contigs} nodes.")
    table = BBHashTable()
    table.initialize(all_kmers)

    # build tables linking:
    # * mphf hash to k-mer hash (for checking exactness)
    # * mphf hash to cDBG ID
    # * cDBG ID to node size (in k-mers)

    print("second pass.")
    records_iter = records_iter_fn()
    sizes = {}
    for n, record in enumerate(records_iter):
        if n % 50000 == 0 and n:
            print("... contig {} of {}".format(n, n_contigs), end="\r")

        # node ID is record name, must go from 0 to total-1
        cdbg_id = int(record.name)
        # get 64-bit numbers for each k-mer
        kmers = hash_sequence(record.sequence, ksize)

        # for each k-mer, find its MPHF hashval, & link to info.
        for kmer in kmers:
            table[kmer] = cdbg_id

        sizes[cdbg_id] = len(kmers)

    print(f"loaded {n_contigs} contigs in pass2.\n")
    assert n == max(table.mphf_to_value), (n, max(table.mphf_to_value))

    return table, sizes


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument("catlas_prefix")
    p.add_argument("--contigs-db", required=True)
    p.add_argument("-k", "--ksize", default=31, type=int)
    a = p.parse_args(argv)

    mphf_filename = os.path.join(a.catlas_prefix, "contigs.mphf")
    array_filename = os.path.join(a.catlas_prefix, "contigs.indices")
    sizes_filename = os.path.join(a.catlas_prefix, "contigs.sizes")

    sqlite_db = sqlite3.connect(a.contigs_db)

    def create_records_iter():
        from spacegraphcats.search import search_utils

        print(f"reading cDBG nodes from sqlite DB {a.contigs_db}")
        return search_utils.contigs_iter_sqlite(sqlite_db)

    table, sizes = build_mphf(a.ksize, create_records_iter)
    print(f"done! saving to {mphf_filename} and {array_filename}")

    table.save(mphf_filename, array_filename)
    with open(sizes_filename, "wb") as fp:
        pickle.dump(sizes, fp)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
