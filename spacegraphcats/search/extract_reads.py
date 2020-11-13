#! /usr/bin/env python
"""
Retrieve the reads for a list of cDBG nodes.  Consumes the output of
extract_nodes_by_query to get the list of nodes, and then uses the
labeled cDBG output by .cdbg.label_cdbg to find reads that overlap with
the unitigs in those nodes.
"""
import argparse
import os
import sys
import time
import gzip

from .search_utils import get_reads_by_cdbg


def main(argv=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument("readsfile")
    p.add_argument("labeled_reads_sqlite")
    p.add_argument("node_list_file", help="a cdbg_ids.txt.gz file")
    p.add_argument("-o", "--output")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    dbfilename = args.labeled_reads_sqlite
    if not os.path.exists(dbfilename):
        print("sqlite file {} does not exist".format(dbfilename))
        sys.exit(-1)

    if not args.output:
        outname = args.node_list_file + ".reads.gz"
    else:
        outname = args.output

    if outname.endswith(".gz"):
        outfp = gzip.open(outname, "wt")
    else:
        outfp = open(outname, "wt")

    with gzip.open(args.node_list_file, "rt") as fp:
        cdbg_shadow = set([int(x.strip()) for x in fp if x.strip()])

    print("extracting reads to {}.".format(outname))

    start = time.time()

    total_bp = 0
    total_seqs = 0
    is_fastq = 0

    print("loading labels and querying for matching sequences..")
    reads_iter = get_reads_by_cdbg(dbfilename, args.readsfile, cdbg_shadow)
    for n, (record, offset_f) in enumerate(reads_iter):
        if n % 10000 == 0:
            print("...at n {} ({:.1f}% of file)".format(n, offset_f * 100), end="\r")

            if n == 0:
                if hasattr(record, "quality"):
                    is_fastq = 1

        # output!
        if is_fastq:
            outfp.write(
                "@{}\n{}\n+\n{}\n".format(record.name, record.sequence, record.quality)
            )
        else:
            outfp.write(">{}\n{}\n".format(record.name, record.sequence))
        total_bp += len(record.sequence)
        total_seqs += 1

    print("")
    print("fetched {} reads, {} bp matching nodes.".format(total_seqs, total_bp))
    end = time.time()
    print(
        "total read retrieval time (including database query): {:.2f}s".format(
            end - start
        )
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
