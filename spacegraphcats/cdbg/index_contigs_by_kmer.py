"""
Use Minimal Perfect Hashing (see BBHash) to construct a fast lookup
table connecting k-mers in the cDBG to cDBG node IDs.

Input: a directory containing a contigs.fa.gz

Output: contigs.fa.gz.mphf, a BBHash MPHF savefile; and contigs.fa.gz.indices,
a numpy savez file containing mphf_to_kmer, kmer_to_cdbg, and sizes.

Note: relies on the fact that for a cDBG constructed at a particular k,
no k-mer will appear in more than one cDBG node and every k-mer will
be in at least one cDBG node, i.e. k <-> cdbg_id is bijective.
"""
import sys
import os
import screed
import khmer
import argparse
import bbhash
import numpy


hashing_fn = None
hashing_ksize = None


def hash_sequence(seq, ksize):
    global hashing_fn, hashing_ksize
    if hashing_fn is None:
        kh = khmer.Nodetable(ksize, 1, 1)
        hashing_fn = kh.get_kmer_hashes
        hashing_ksize = ksize

    assert hashing_ksize == ksize
    return hashing_fn(seq)


class MPHF_KmerIndex(object):
    """
    Support kmer -> cDBG node id queries.
    """

    def __init__(self, mphf, mphf_to_kmer, mphf_to_cdbg_id, cdbg_sizes):
        self.mphf = mphf
        self.mphf_to_kmer = mphf_to_kmer
        self.mphf_to_cdbg_id = mphf_to_cdbg_id
        self.cdbg_sizes = cdbg_sizes

    def get_cdbg_size(self, cdbg_id):
        "Return the size (in k-mers) of the given cDBG ID."
        return self.cdbg_sizes[int(cdbg_id)]

    def get_cdbg_id(self, kmer_hash):
        # do we have this kmer in the mphf?
        mphf_hash = self.mphf.lookup(kmer_hash)
        if mphf_hash:
            # double check that mphf maps to correct hash
            if self.mphf_to_kmer[mphf_hash] == kmer_hash:
                cdbg_id = self.mphf_to_cdbg_id[mphf_hash]
                return cdbg_id
        return None

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

    def count_cdbg_matches(self, query_kmers, verbose=True):
        "Return a dictionary containing cdbg_id -> # of matches in query_kmers"
        cdbg_matches = {}
        n_matched = 0
        n = 0
        for n, hashval in enumerate(query_kmers):
            if n and n % 1000000 == 0 and verbose:
                print("matching ...", n, end="\r")

            cdbg_id = self.get_cdbg_id(hashval)
            if cdbg_id is not None:
                cdbg_matches[cdbg_id] = cdbg_matches.get(cdbg_id, 0) + 1
                n_matched += 1

        if n and verbose:
            print("... found {} matches to {} k-mers total.".format(n_matched, n + 1))

        return cdbg_matches

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
        mphf_filename = os.path.join(catlas_prefix, "contigs.fa.gz.mphf")
        array_filename = os.path.join(catlas_prefix, "contigs.fa.gz.indices")
        if not os.path.exists(mphf_filename):
            raise FileNotFoundError(mphf_filename)
        mphf = bbhash.load_mphf(mphf_filename)
        with open(array_filename, "rb") as fp:
            np_dict = numpy.load(fp)
            mphf_to_kmer = np_dict["mphf_to_kmer"]
            mphf_to_cdbg = np_dict["kmer_to_cdbg"]
            cdbg_sizes = np_dict["sizes"]
        return cls(mphf, mphf_to_kmer, mphf_to_cdbg, cdbg_sizes)


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
    print("loaded {} contigs.\n".format(n_contigs))

    # build MPHF (this is the CPU intensive bit)
    print("building MPHF for {} k-mers in {} nodes.".format(len(all_kmers), n_contigs))
    x = bbhash.PyMPHF(all_kmers, len(all_kmers), 4, 1.0)

    # build tables linking:
    # * mphf hash to k-mer hash (for checking exactness)
    # * mphf hash to cDBG ID
    # * cDBG ID to node size (in k-mers)

    mphf_to_kmer = numpy.zeros(len(all_kmers), numpy.uint64)
    mphf_to_cdbg = numpy.zeros(len(all_kmers), numpy.uint32)
    sizes = numpy.zeros(n_contigs, numpy.uint32)

    print("second pass.")
    records_iter = records_iter_fn()
    for n, record in enumerate(records_iter):
        if n % 50000 == 0 and n:
            print("... contig {} of {}".format(n, n_contigs), end="\r")

        # node ID is record name, must go from 0 to total-1
        cdbg_id = int(record.name)

        # get 64-bit numbers for each k-mer (doesn't really matter what hash)
        kmers = hash_sequence(record.sequence, ksize)

        # for each k-mer, find its MPHF hashval, & link to info.
        for kmer in kmers:
            mphf = x.lookup(kmer)
            mphf_to_kmer[mphf] = kmer
            mphf_to_cdbg[mphf] = cdbg_id

        # record each node size, while we're here.
        sizes[cdbg_id] = len(kmers)

    print("loaded {} contigs in pass2.\n".format(n_contigs))
    assert n == max(mphf_to_cdbg), (n, max(mphf_to_cdbg))

    return x, mphf_to_kmer, mphf_to_cdbg, sizes


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument("catlas_prefix")
    p.add_argument("-k", "--ksize", default=31, type=int)
    a = p.parse_args(argv)

    contigs_filename = os.path.join(a.catlas_prefix, "contigs.fa.gz")
    mphf_filename = os.path.join(a.catlas_prefix, "contigs.fa.gz.mphf")
    array_filename = os.path.join(a.catlas_prefix, "contigs.fa.gz.indices")

    def create_records_iter():
        print("reading cDBG nodes from {}".format(contigs_filename))
        return screed.open(contigs_filename)

    x, mphf_to_kmer, mphf_to_cdbg, sizes = build_mphf(a.ksize, create_records_iter)

    print("done! saving to {} and {}".format(mphf_filename, array_filename))

    x.save(mphf_filename)
    with open(array_filename, "wb") as fp:
        numpy.savez_compressed(
            fp, mphf_to_kmer=mphf_to_kmer, kmer_to_cdbg=mphf_to_cdbg, sizes=sizes
        )


if __name__ == "__main__":
    main(sys.argv[1:])
