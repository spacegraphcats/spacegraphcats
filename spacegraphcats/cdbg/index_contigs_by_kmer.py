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


def build_mphf(kh, records_iter_fn):
    # build a list of all k-mers in the cDBG
    all_kmers = list()

    records_iter = records_iter_fn()
    for n, record in enumerate(records_iter):
        if n % 50000 == 0 and n:
            print('... contig', n, end='\r')

        kmers = kh.get_kmer_hashes(record.sequence)
        all_kmers.extend(list(kmers))

    n_contigs = n + 1
    print('loaded {} contigs.\n'.format(n_contigs))

    # build MPHF (this is the CPU intensive bit)
    print('building MPHF for {} k-mers in {} nodes.'.format(len(all_kmers), n_contigs))
    x = bbhash.PyMPHF(all_kmers, len(all_kmers), 4, 1.0)

    # build tables linking:
    # * mphf hash to k-mer hash (for checking exactness)
    # * mphf hash to cDBG ID
    # * cDBG ID to node size (in k-mers)

    mphf_to_kmer = numpy.zeros(len(all_kmers), numpy.uint64)
    mphf_to_cdbg = numpy.zeros(len(all_kmers), numpy.uint32)
    sizes = numpy.zeros(n_contigs, numpy.uint32)

    print('second pass.')
    records_iter = records_iter_fn()
    for n, record in enumerate(records_iter):
        if n % 50000 == 0 and n:
            print('... contig {} of {}'.format(n, n_contigs), end='\r')

        # node ID is record name, must go from 0 to total-1
        cdbg_id = int(record.name)

        # get 64-bit numbers for each k-mer (doesn't really matter what hash)
        kmers = kh.get_kmer_hashes(record.sequence)

        # for each k-mer, find its MPHF hashval, & link to info.
        for kmer in kmers:
            mphf = x.lookup(kmer)
            mphf_to_kmer[mphf] = kmer
            mphf_to_cdbg[mphf] = cdbg_id

        # record each node size, while we're here.
        sizes[cdbg_id] = len(kmers)

    print('loaded {} contigs in pass2.\n'.format(n_contigs))
    assert n == max(mphf_to_cdbg), (n, max(mphf_to_cdbg))

    return x, mphf_to_kmer, mphf_to_cdbg, sizes


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix')
    p.add_argument('-k', '--ksize', default=31, type=int)
    a = p.parse_args(argv)

    kh = khmer.Nodetable(a.ksize, 1, 1)

    contigs_filename = os.path.join(a.catlas_prefix, 'contigs.fa.gz')
    mphf_filename = os.path.join(a.catlas_prefix, 'contigs.fa.gz.mphf')
    array_filename = os.path.join(a.catlas_prefix, 'contigs.fa.gz.indices')

    def create_records_iter():
        print('reading cDBG nodes from {}'.format(contigs_filename))
        return screed.open(contigs_filename)

    x, mphf_to_kmer, mphf_to_cdbg, sizes = build_mphf(kh, create_records_iter)

    print('done! saving to {} and {}'.format(mphf_filename, array_filename))

    x.save(mphf_filename)
    with open(array_filename, 'wb') as fp:
        numpy.savez_compressed(fp,
                               mphf_to_kmer=mphf_to_kmer,
                               kmer_to_cdbg=mphf_to_cdbg,
                               sizes=sizes)


if __name__ == '__main__':
    main(sys.argv[1:])
