#! /usr/bin/env python
"""
Given a catlas and a sample, calculate the abundances of the catlas domset
in the sample.

One use for this is when you have a catlas built from multiple samples,
and you want to output expression directly on the domset using
abundances from individual samples.
"""
import sys
import argparse
import time
from collections import defaultdict
import csv

import spacegraphcats
from spacegraphcats.utils.logging import notify, error, debug
from spacegraphcats.search import search_utils
from spacegraphcats.search import MPHF_KmerIndex, hash_sequence
from spacegraphcats.search.catlas import CAtlas
import screed


def main():
    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument('cdbg_prefix', help='cdbg prefix')
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('samples', nargs='+')
    p.add_argument('-k', '--ksize', default=31, type=int,
                   help='k-mer size (default: 31)')
    args = p.parse_args()

    # load catlas DAG
    catlas = CAtlas(args.cdbg_prefix, args.catlas_prefix, load_sizefile=True)
    notify('loaded {} nodes from catlas {}', len(catlas), args.catlas_prefix)
    notify('loaded {} layer 1 catlas nodes', len(catlas.layer1_to_cdbg))

    ki_start = time.time()
    kmer_idx = MPHF_KmerIndex.from_directory(args.cdbg_prefix)
    notify('loaded {} k-mers in index ({:.1f}s)',
           len(kmer_idx), time.time() - ki_start)

    # calculate the k-mer sizes beneath node, for each catlas node.
    catlas.decorate_with_index_sizes(kmer_idx)

    #
    # run through each sample, calculating all the things.
    #

    for sample in args.samples:
        # here, cdbg_counts is the summed k-mer abundance from the counts
        # of k-mers in the sample, for each cdbg_id.

        notify(f'reading sequences from {sample}')
        cdbg_counts = defaultdict(int)
        total_hashes = 0
        for record in screed.open(sample):
            hashes = hash_sequence(record.sequence, args.ksize)
            total_hashes += len(hashes)
            for hashval in hashes:
                cdbg_id = kmer_idx.get_cdbg_id(hashval)
                cdbg_counts[cdbg_id] += 1

        notify(f'done! read {total_hashes} total k-mers/abundances.')

        # aggregate to layer1 domset across all cDBG nodes
        dom_counts = defaultdict(int)
        dom_sizes = defaultdict(int)
        for cdbg_id, count in catlas.cdbg_sizes.items():
            count = cdbg_counts[cdbg_id]
            dom_id = catlas.cdbg_to_layer1[cdbg_id]
            dom_counts[dom_id] += count
            dom_sizes[dom_id] += catlas.cdbg_sizes[cdbg_id]

        outfile = f"{sample}.cdbg_abund.csv"
        notify(f'outputting cDBG abundances to {outfile}')
        total_cdbg_counts = 0
        with open(outfile, 'wt') as fp:
            w = csv.writer(fp)
            w.writerow(['cdbg_id', 'abund', 'node_size'])
            for cdbg_id in sorted(catlas.cdbg_sizes):   # get all cDBG IDs
                count = cdbg_counts[cdbg_id]
                size = catlas.cdbg_sizes[cdbg_id]
                w.writerow([cdbg_id, count, size])
                total_cdbg_counts += count

        # track counts per level, for double-checking
        level_counts = defaultdict(int)
        level_sizes = defaultdict(int)

        outfile = f"{sample}.dom_abund.csv"
        notify(f'outputting dom node abundances for all levels to {outfile}')
        with open(outfile, 'wt') as fp:
            w = csv.writer(fp)
            w.writerow(['dom_id', 'abund', 'size'])
            for node_id in sorted(catlas):
                leaves = catlas.leaves([node_id])

                count = 0
                size = 0
                for leaf_id in leaves:
                    count += dom_counts[leaf_id]
                    size += dom_sizes[leaf_id]
                w.writerow([node_id, count, size])

                level = catlas.levels[node_id]
                level_counts[level] += count
                level_sizes[level] += size

        print('level_counts:', level_counts)
        print('total_cdbg_counts:', total_cdbg_counts)
        print('dom_counts:', sum(dom_counts.values()))

        assert total_cdbg_counts == sum(dom_counts.values())

        print('total_cdbg_sizes:', sum(catlas.cdbg_sizes.values()))
        print('sum level sizes:', level_sizes)


    return 0


if __name__ == '__main__':
    sys.exit(main())
