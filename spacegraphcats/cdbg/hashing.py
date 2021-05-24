"""
Hashing code for k-mers -> hash values.
"""
import os
import pickle

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


class MPHF_KmerIndex:
    """
    Support kmer -> cDBG node id queries.
    """

    def __init__(self, ksize, bbhash_table, sizes):
        self.ksize = ksize
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
                # propagate up the catlas.
                sub_size = 0
                for child_id in catlas.children[node_id]:
                    sub_size += dom_matches.get(child_id, 0)

                if sub_size:
                    dom_matches[node_id] = sub_size

        return dom_matches

    def save(self, location):
        "Save kmer index."
        mphf_filename = os.path.join(location, "contigs.mphf")
        array_filename = os.path.join(location, "contigs.indices")
        sizes_filename = os.path.join(location, "contigs.sizes")

        self.table.save(mphf_filename, array_filename)
        with open(sizes_filename, "wb") as fp:
            pickle.dump((self.ksize, self.sizes), fp)

    @classmethod
    def from_directory(cls, location):
        "Load kmer index created by search.contigs_by_kmer."
        mphf_filename = os.path.join(location, "contigs.mphf")
        array_filename = os.path.join(location, "contigs.indices")
        sizes_filename = os.path.join(location, "contigs.sizes")

        table = BBHashTable.load(mphf_filename, array_filename)
        with open(sizes_filename, "rb") as fp:
            ksize, sizes = pickle.load(fp)

        return cls(ksize, table, sizes)
