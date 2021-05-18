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
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('samples', nargs='+')
    p.add_argument('-k', '--ksize', default=31, type=int,
                   help='k-mer size (default: 31)')
    args = p.parse_args()

    # load catlas DAG
    catlas = CAtlas(args.catlas_prefix, load_sizefile=True)
    notify('loaded {} nodes from catlas {}', len(catlas), args.catlas_prefix)
    notify('loaded {} layer 1 catlas nodes', len(catlas.layer1_to_cdbg))

    ki_start = time.time()
    kmer_idx = MPHF_KmerIndex.from_catlas_directory(args.catlas_prefix)
    notify('loaded {} k-mers in index ({:.1f}s)',
           len(kmer_idx), time.time() - ki_start)

    # calculate the k-mer sizes for each catlas node.
    catlas.decorate_with_index_sizes(kmer_idx)

    for sample in args.samples:
        notify(f'reading sequences from {sample}')
        cdbg_counts = defaultdict(int)
        total_hashes = 0
        for record in screed.open(sample):
            hashes = hash_sequence(record.sequence, args.ksize)
            total_hashes += len(hashes)
            for hashval in hashes:
                cdbg_id = kmer_idx.get_cdbg_id(hashval)
                cdbg_counts[cdbg_id] += 1

        notify(f'done! read {total_hashes} total k-mers.')

        # aggregate to domset
        dom_counts = defaultdict(int)
        for cdbg_id, count in cdbg_counts.items():
            if cdbg_id is None:
                continue
            dom_id = catlas.cdbg_to_layer1[cdbg_id]
            dom_counts[dom_id] += count

        outfile = f"{sample}.cdbg_abund.csv"
        notify(f'outputting cDBG abundances to {outfile}')
        with open(outfile, 'wt') as fp:
            w = csv.writer(fp)
            w.writerow(['cdbg_id', 'abund'])
            for cdbg_id in sorted(catlas.cdbg_sizes):   # get all cDBG IDs
                count = cdbg_counts[cdbg_id]
                w.writerow([cdbg_id, count])

        outfile = f"{sample}.dom_abund.csv"
        notify(f'outputting dom node abundances to {outfile}')
        with open(outfile, 'wt') as fp:
            w = csv.writer(fp)
            w.writerow(['dom_id', 'abund'])
            for node_id in sorted(catlas):
                if catlas.levels[node_id] == 1:
                    count = dom_counts[node_id]
                    w.writerow([node_id, count])

    return 0


if __name__ == '__main__':
    sys.exit(main())
