#! /usr/bin/env python
"""
A simple DNA-kmer-based min-set-cov implementation: find the minimum
set cover from reads for a reference set of k-mers.

The greedy implementation used here will prioritize longer reads over shorter
reads.

Based on sourmash gather.
"""
import sys
import argparse

import screed
import sourmash
from spacegraphcats.utils.counter_gather import CounterGather


def main(argv):
    """\
    Use k-mers to construct a minimum set cover for (small) reference sets
    using (somewhat larger?) collections of reads.

    Note, everything is kept in memory here, so apply ...wisely.
    """
    p = argparse.ArgumentParser(description=main.__doc__)

    p.add_argument("reference", help="reference contigs or reads")
    p.add_argument("reads", help="read file(s)", nargs='+')
    p.add_argument("-o", "--output", required=True, help="output file")

    p.add_argument("-k", "--ksize", default=21, type=int,
        help="DNA k-mer size (default: 10)"
    )

    args = p.parse_args()

    ksize = args.ksize
    mh = sourmash.MinHash(n=0, scaled=1, ksize=ksize)
    print(f"Using DNA k-mer size {ksize}")

    print(f"Loading reference k-mers from '{args.reference}'")
    ref_kmers = set()
    for record in screed.open(args.reference):
        kmers = mh.seq_to_hashes(record.sequence, force=True)
        ref_kmers.update(kmers)

    print(f"Loaded {len(ref_kmers)} reference k-mers total.")
    counter = CounterGather(ref_kmers)

    found_kmers = set()
    reads_dict = {}
    n = 0
    for filename in args.reads:
        print(f"Loading reads from '{filename}'")
        for record in screed.open(filename):
            kmers = set(mh.seq_to_hashes(record.sequence, force=True))
            intersect_kmers = kmers & ref_kmers
            if intersect_kmers:
                reads_dict[n] = record
                counter.add(kmers, n)
                found_kmers.update(intersect_kmers)
                n += 1

    print(f"Found {len(reads_dict)} matching reads, total.")
    pcnt = len(found_kmers) / len(ref_kmers) * 100
    print(f"{len(found_kmers)} kmers found of {len(ref_kmers)} total ({pcnt:.1f}%)")
    print("Running gather...")

    remaining_kmers = set(ref_kmers)
    filtered_names = {}
    x = counter.peek(remaining_kmers)
    while x:
        cont, match_set, name, intersect_set = x
        counter.consume(intersect_set)
        remaining_kmers -= intersect_set
        assert name not in filtered_names
        filtered_names[name] = len(intersect_set)

        x = counter.peek(remaining_kmers)

    print(f"...minimum set cover is {len(filtered_names)} reads.")

    total_bp = 0
    with open(args.output, 'wt') as fp:
        for read_id in filtered_names:
            read = reads_dict[read_id]
            fp.write(f'>{read.name}\n{read.sequence}\n')
            total_bp += len(record.sequence)

    print(f"output {total_bp} bp to '{args.output}'")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
