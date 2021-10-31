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
import os

from spacegraphcats.utils.logging import notify
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
    p.add_argument('--outdir', default=None,
                   help='output directory; by default, {catlas_prefix}_abund')
    p.add_argument('--prefix', default="",
                   help="A prefix to use for all output file names.")
    args = p.parse_args()

    # load catlas DAG
    catlas = CAtlas(args.cdbg_prefix, args.catlas_prefix, load_sizefile=True)
    notify('loaded {} nodes from catlas {}', len(catlas), args.catlas_prefix)
    notify('loaded {} layer 1 catlas nodes', len(catlas.layer1_to_cdbg))

    ki_start = time.time()
    kmer_idx = MPHF_KmerIndex.from_directory(args.cdbg_prefix)
    notify('loaded {} k-mers in index ({:.1f}s)',
           len(kmer_idx), time.time() - ki_start)

    catlas_root_children = set(catlas.children[catlas.root])
    max_level = max(catlas.levels.values())

    # calculate the k-mer sizes beneath node, for each catlas node.
    catlas.decorate_with_index_sizes(kmer_idx)

    # set and/or create output dir
    outdir = args.outdir
    if not outdir:
        catlas_name = os.path.dirname(args.catlas_prefix)
        outdir = f"{catlas_name}_abund"
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    notify(f"Outputting files to directory '{outdir}'")

    #
    # run through each sample, calculating all the things.
    #

    for sample in args.samples:
        sample_name = os.path.basename(sample)

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
        notify('--')

        # aggregate to layer1 domset across all cDBG nodes
        dom_counts = defaultdict(int)
        dom_sizes = defaultdict(int)
        for cdbg_id, count in catlas.cdbg_sizes.items():
            count = cdbg_counts[cdbg_id]
            dom_id = catlas.cdbg_to_layer1[cdbg_id]
            dom_counts[dom_id] += count
            dom_sizes[dom_id] += catlas.cdbg_sizes[cdbg_id]

        outfile = f"{args.prefix}{sample_name}.cdbg_abund.csv"
        outfile = os.path.join(outdir, outfile)

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

        outfile = f"{args.prefix}{sample_name}.dom_abund.csv"
        outfile = os.path.join(outdir, outfile)

        notify(f'outputting dom node abundances for all levels to {outfile}')
        with open(outfile, 'wt') as fp:
            w = csv.writer(fp)
            w.writerow(['dom_id', 'level', 'abund', 'size'])

            for node_id in sorted(catlas):
                leaves = catlas.leaves([node_id])
                level = catlas.levels[node_id]

                count = 0
                size = 0
                for leaf_id in leaves:
                    count += dom_counts[leaf_id]
                    size += dom_sizes[leaf_id]
                w.writerow([node_id, level, count, size])

                level_counts[level] += count
                level_sizes[level] += size

            # recover level info for those nodes that represent entire
            # components of the graph, and hence are attached
            # directly to catlas root; they no longer appear in catlas
            # structure (until root node is reached).
            #
            # NOTE: these node_ids will appear multiple times in the
            # dom_abund.csv file!

            for node_id in catlas_root_children:
                stop_level = catlas.levels[node_id]
                leaves = catlas.leaves([node_id])

                count = 0
                size = 0
                for leaf_id in leaves:
                    count += dom_counts[leaf_id]
                    size += dom_sizes[leaf_id]

                for level in range(stop_level + 1, max_level):
                    level_counts[level] += count
                    level_sizes[level] += size
                    w.writerow([node_id, level, count, size])

        level_count_values = set(level_counts.values())
        assert len(level_count_values) == 1
        level_count = list(level_count_values)[0]

        level_sizes_values = set(level_sizes.values())
        assert len(level_sizes_values) == 1
        level_size = list(level_sizes_values)[0]

        print(f'- each catlas level has {level_count} sum counts')
        print('- total_cdbg_counts:', total_cdbg_counts)
        print('- dom_counts:', sum(dom_counts.values()))

        assert total_cdbg_counts == sum(dom_counts.values())

        print(f'- each catlas level has {level_size} k-mers underneath')
        print('- total_cdbg_sizes:', sum(catlas.cdbg_sizes.values()))

    return 0


if __name__ == '__main__':
    sys.exit(main())
