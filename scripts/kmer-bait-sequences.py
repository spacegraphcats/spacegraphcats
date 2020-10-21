#! /usr/bin/env python
"""
Use a set of query files to extract sequences with k-mer overlaps.
"""
import sys
import screed
import argparse
import khmer

REPORTING_BP = 1e6


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("extract_from")
    parser.add_argument("queries", nargs="+")
    parser.add_argument("-o", "--output")
    parser.add_argument("-k", "--ksize", default=31, type=int)
    args = parser.parse_args()

    assert args.output, "must provide an argument to -o"

    kh = khmer.Nodetable(args.ksize, 1, 1)

    query_kmers = set()
    for queryfile in args.queries:
        print("loading query", queryfile, file=sys.stderr)
        for record in screed.open(queryfile):
            query_kmers.update(kh.get_kmer_hashes(record.sequence))

    print("loaded {} k-mers".format(len(query_kmers)), file=sys.stderr)

    bp = 0
    threshold = REPORTING_BP
    m = 0
    n = 0
    bp_out = 0

    with open(args.output, "wt") as fp:
        print("searching sequences in {}".format(args.extract_from))
        for record in screed.open(args.extract_from):
            n += 1
            bp += len(record.sequence)
            if bp > threshold:
                threshold += REPORTING_BP
                print("... read {} bp in {} sequences".format(bp, n), file=sys.stderr)
                print(
                    "==> found {} bp in {} sequences".format(bp_out, m), file=sys.stderr
                )

            x = set(kh.get_kmer_hashes(record.sequence))
            if x.intersection(query_kmers):
                fp.write(">{}\n{}\n".format(record.name, record.sequence))
                bp_out += len(record.sequence)
                m += 1

    print("... read {} bp in {} sequences".format(bp, n), file=sys.stderr)
    print("==> found {} bp in {} sequences".format(bp_out, m), file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
