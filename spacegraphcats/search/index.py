import os
import bbhash
import numpy


class MPHF_KmerIndex(object):
    """
    Support kmer -> cDBG node id queries enabled by search.contigs_by_kmer.
    """
    def __init__(self, mphf, mphf_to_kmer, mphf_to_cdbg_id, cdbg_sizes):
        self.mphf = mphf
        self.mphf_to_kmer = mphf_to_kmer
        self.mphf_to_cdbg_id = mphf_to_cdbg_id
        self.cdbg_sizes = cdbg_sizes

    def get_cdbg_size(self, cdbg_id):
        return self.cdbg_sizes[cdbg_id]

    def get_cdbg_id(self, kmer_hash):
        mphf_hash = self.mphf.lookup(kmer_hash)
        if mphf_hash is None:
            return None
        elif self.mphf_to_kmer[mphf_hash] == kmer_hash:
            cdbg_id = self.mphf_to_cdbg_id[mphf_hash]
            return cdbg_id
        return None

    def build_catlas_node_sizes(self, catlas):
        """Count the number of kmers in the shadow of each catlas node."""
        node_kmer_sizes = {}
        for node_id in catlas:
            level = catlas.levels[node_id]
            if level == 1:
                total_kmers = 0
                for cdbg_node in catlas.layer1_to_cdbg.get(node_id):
                    total_kmers += self.get_cdbg_size(cdbg_node)

                node_kmer_sizes[node_id] = total_kmers
            else:
                sub_size = 0
                for child_id in catlas.children[node_id]:
                    sub_size += node_kmer_sizes[child_id]
                node_kmer_sizes[node_id] = sub_size

        return node_kmer_sizes

    def get_match_counts(self, query_kmers):
        "Return a dictionary containing cdbg_id -> # of matches in query_kmers"
        match_counts = {}
        n_matched = 0
        n = 0
        for n, hashval in enumerate(query_kmers):
            if n % 1000000 == 0:
                print('matching ...', n, end='\r')

            cdbg_id = self.get_cdbg_id(hashval)
            if cdbg_id is not None:
                match_counts[cdbg_id] = match_counts.get(cdbg_id, 0) + 1
                n_matched += 1

        if n:
            print('... found {} matches to {} k-mers total.'.format(n_matched,
                                                                    n + 1))

        return match_counts

    def build_catlas_match_counts(self, match_counts, catlas):
        node_match_counts = {}
        total_kmers_in_query_nodes = 0.
        for node_id in catlas:
            level = catlas.levels[node_id]
            if level == 1:
                query_kmers = 0
                all_kmers_in_node = 0
                for cdbg_node in catlas.layer1_to_cdbg.get(node_id):
                    query_kmers += match_counts.get(cdbg_node, 0)
                    all_kmers_in_node += self.get_cdbg_size(cdbg_node)

                if query_kmers:
                    node_match_counts[node_id] = query_kmers
                    total_kmers_in_query_nodes += all_kmers_in_node

                    assert all_kmers_in_node >= query_kmers

            else:
                sub_size = 0
                for child_id in catlas.children[node_id]:
                    sub_size += node_match_counts.get(child_id, 0)

                if sub_size:
                    node_match_counts[node_id] = sub_size

        return node_match_counts

    @classmethod
    def from_catlas_directory(cls, catlas_prefix):
        "Load kmer index created by search.contigs_by_kmer."
        mphf_filename = os.path.join(catlas_prefix, 'contigs.fa.gz.mphf')
        array_filename = os.path.join(catlas_prefix, 'contigs.fa.gz.indices')
        mphf = bbhash.load_mphf(mphf_filename)
        with open(array_filename, 'rb') as fp:
            np_dict = numpy.load(fp)
            mphf_to_kmer = np_dict['mphf_to_kmer']
            mphf_to_cdbg = np_dict['kmer_to_cdbg']
            cdbg_sizes = np_dict['sizes']
        return cls(mphf, mphf_to_kmer, mphf_to_cdbg, cdbg_sizes)
