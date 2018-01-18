import sys
import os
import screed
import khmer
import argparse
import bbhash
import numpy


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix')
    a = p.parse_args()
    
    kh = khmer.Nodetable(31, 1, 1)

    contigs_filename = os.path.join(a.catlas_prefix, 'contigs.fa.gz')
    mphf_filename = os.path.join(a.catlas_prefix, 'contigs.fa.gz.mphf')
    array_filename = os.path.join(a.catlas_prefix, 'contigs.fa.gz.indices')

    d = {}
    sizes = {}

    all_kmers = list()
    print('reading cDBG nodes from {}'.format(contigs_filename))
    for n, record in enumerate(screed.open(contigs_filename)):
        if n % 50000 == 0 and n:
            print('... contig', n)

        kmers = kh.get_kmer_hashes(record.sequence)
        all_kmers.extend(list(kmers))

    n_contigs = n + 1

    x = bbhash.PyMPHF(all_kmers, len(all_kmers), 4, 1.0)

    mphf_to_kmer = numpy.zeros(len(all_kmers), numpy.uint64)
    kmer_to_cdbg = numpy.zeros(len(all_kmers), numpy.uint32)
    sizes = numpy.zeros(n_contigs, numpy.uint32)

    print('second pass; reading cDBG nodes from {}'.format(contigs_filename))
    for n, record in enumerate(screed.open(contigs_filename)):
        if n % 50000 == 0 and n:
            print('... contig {} of {}'.format(n, n_contigs))

        kmers = kh.get_kmer_hashes(record.sequence)
        for kmer in kmers:
            ph = x.lookup(kmer)
            mphf_to_kmer[ph] = kmer
            kmer_to_cdbg[ph] = int(record.name)

        sizes[n] = len(kmers)

    x.save(mphf_filename)
    with open(array_filename, 'wb') as fp:
        numpy.savez_compressed(fp,
                               mphf_to_kmer=mphf_to_kmer,
                               kmer_to_cdbg=kmer_to_cdbg,
                               sizes=sizes)


if __name__ == '__main__':
    main()
