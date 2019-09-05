#! /usr/bin/env python
"""
Retrieve nodes by MinHash hashval, using the indices created by
spacegraphcats.cdbg.index_cdbg_by_minhash.
"""
import argparse
import csv
import gzip
import os
import sys
import time
import pickle

import khmer
import screed
import sourmash_lib
from sourmash_lib import MinHash

from . import search_utils
from .index import MPHF_KmerIndex
from .catlas import CAtlas


class QueryOutput:
    def __init__(self, query_hashval, catlas, cdbg_shadow):
        self.query_hashval = query_hashval
        self.catlas = catlas
        self.cdbg_shadow = cdbg_shadow
        self.total_bp = 0
        self.total_seq = 0
        self.contigs = []

    def __add_sequence(self, sequence):
        self.total_bp += len(sequence)
        self.total_seq += 1

    def retrieve_contigs(self, contigs):
        "extract contigs using cDBG shadow."

        # node list.
        # track extracted info
        retrieve_start = time.time()
        # walk through the contigs, retrieving.
        print('extracting contigs...')
        for n, record in enumerate(
                search_utils.get_contigs_by_cdbg(contigs,
                                                 self.cdbg_shadow)):
            self.total_seq += 1
            self.total_bp += len(record.sequence)

            if n and n % 10000 == 0:
                offset_f = self.total_seq / len(self.cdbg_shadow)
                print('...at n {} ({:.1f}% of shadow)'.format(self.total_seq,
                      offset_f * 100), end='\r')

            self.contigs.append((record.name, record.sequence))

        # done - got all contigs!
        print('...fetched {} contigs, {} bp matching combined frontiers. '
              ' ({:.1f}s)'.format(self.total_seq, self.total_bp,
                                  time.time() - retrieve_start))


    def write(self, csv_writer, csvoutfp, outdir):
        hashval = self.query_hashval
        bp = self.total_bp
        seqs = self.total_seq
        k = 0

        # output to results.csv!
        csv_writer.writerow([hashval, bp, seqs, k])
        csvoutfp.flush()

        # write out cDBG IDs
        q_name = str(hashval)
        cdbg_listname = os.path.basename(q_name) + '.cdbg_ids.txt.gz'
        with gzip.open(os.path.join(outdir, cdbg_listname), 'wt') as fp:
            fp.write("\n".join([str(x) for x in sorted(self.cdbg_shadow)]))

        # write out contigs
        contigs_outname = os.path.basename(q_name) + '.contigs.fa.gz'
        with gzip.open(os.path.join(outdir, contigs_outname), 'wt') as fp:
            for name, sequence in self.contigs:
                fp.write('>{}\n{}\n'.format(name, sequence))


def execute_query(hashval, catlas, hashval_to_contig_id):
    cdbg_node = hashval_to_contig_id.get(hashval, None)

    if not cdbg_node:
        return None

    catlas_node_id = catlas.cdbg_to_layer1[cdbg_node]

    # calculate level 1 nodes for this frontier in the catlas
    leaves = catlas.leaves([catlas_node_id])

    # calculate associated cDBG nodes
    cdbg_shadow = catlas.shadow(leaves)

    print('done searching!')
    return QueryOutput(hashval, catlas, cdbg_shadow)


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('mh_index_picklefile', help='pickled hashval index')
    p.add_argument('hashval_list', help='file with list of hashvals')
    p.add_argument('output')
    p.add_argument('-k', '--ksize', default=31, type=int,
                   help='k-mer size (default: 31)')
    p.add_argument('--scaled', default=1000, type=float,
                   help="scaled value for contigs minhash output")
    p.add_argument('-v', '--verbose', action='store_true')

    args = p.parse_args(argv)

    # create output directory if it doesn't exist.
    outdir = args.output
    try:
        os.mkdir(outdir)
    except OSError:
        pass
    if not os.path.isdir(outdir):
        print('output {} is not a directory'.format(outdir))
        sys.exit(-1)

    # load picklefile
    with open(args.mh_index_picklefile, 'rb') as fp:
        hashval_to_contig_id = pickle.load(fp)

    # load list of desired hashvals
    hashvals = [ int(x.strip()) for x in open(args.hashval_list, 'rt') ]
    hashvals = set(hashvals)

    if not len(hashvals):
        print('No hash values to search!', file=sys.stderr)
        sys.exit(-1)

    # load catlas DAG
    catlas = CAtlas(args.catlas_prefix)
    print('loaded {} nodes from catlas {}'.format(len(catlas), args.catlas_prefix))
    print('loaded {} layer 1 catlas nodes'.format(len(catlas.layer1_to_cdbg)))

    # find the contigs filename
    contigs_file = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

    # get a single ksize & scaled
    ksize = int(args.ksize)
    scaled = int(args.scaled)

    # record command line
    with open(os.path.join(outdir, 'command.txt'), 'wt') as fp:
        fp.write(str(sys.argv))
        fp.write("\n")

    # output results.csv in the output directory:
    csvoutfp = open(os.path.join(outdir, 'hashval_results.csv'), 'wt')
    csv_writer = csv.writer(csvoutfp)
    csv_writer.writerow(['hashval', 'bp', 'contigs', 'ksize'])

    # iterate over each query, do the thing.
    for hashval in hashvals:
        result = execute_query(hashval, catlas, hashval_to_contig_id)
        result.retrieve_contigs(contigs_file)
        result.write(csv_writer, csvoutfp, outdir)
    # @@@
    # end main loop!

    sys.exit(0)


if __name__ == '__main__':
    main(sys.argv[1:])
