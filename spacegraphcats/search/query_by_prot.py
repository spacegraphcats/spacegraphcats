#! /usr/bin/env python
"""Do a search of the unitigs with a protein query.

Usage:

   query-unitigs-prot.py <query_faa_file> <unitigs_sqlite_db> -k <prot_ksize> \
        --out-prefix <prefix>

This will search <unitigs_sqlite_db> with the protein sequence(s) in
<query_faa_file>, and save the resulting cDBG IDs to <prefix>.nodes.gz.

These cdbg IDs can then be used to retrieve neighborhoods with
extract_neighborhoods_by_cdbg_ids.

CTB note: A similar approach can also be used for DNA k-mer search where
k < cDBG ksize.
"""
import sys
import argparse
import os
import sqlite3
import gzip

import screed
import sourmash
from spacegraphcats.search import search_utils


def main(args):
    p = argparse.ArgumentParser()
    p.add_argument("query")
    p.add_argument("unitigs_db")
    p.add_argument("-k", "--ksize", type=int, default=10, help="protein ksize")
    p.add_argument("--out-prefix")
    p.add_argument(
        "--query-is-dna", help="translate query into protein", action="store_true"
    )
    args = p.parse_args(args)

    out_prefix = args.out_prefix
    if not out_prefix:
        out_prefix = os.path.basename(args.query)

    prot_ksize = args.ksize
    dna_ksize = args.ksize * 3

    dna_mh = sourmash.MinHash(n=0, scaled=1000, ksize=dna_ksize)
    prot_mh = sourmash.MinHash(n=0, scaled=100, ksize=prot_ksize, is_protein=True)
    translate_mh = prot_mh.copy_and_clear()

    ## query: for all of the query sequences (which are protein),
    ## translate them into k-mers at protein ksize.
    query_kmers = set()

    if args.query_is_dna:
        add_as_protein = False
    else:
        add_as_protein = True

    for n, record in enumerate(screed.open(args.query)):
        if n and n % 1000 == 0:
            print(f"... {n} {len(query_kmers)} (query)", end="\r", file=sys.stderr)
        these_hashes = translate_mh.seq_to_hashes(
            record.sequence, is_protein=add_as_protein
        )
        query_kmers.update(these_hashes)
    n += 1

    print(
        f"loaded {n} query sequences {len(query_kmers)} prot kmers (query)",
        file=sys.stderr,
    )

    ## now unitigs... for all of the unitigs (which are DNA),
    ## first: translate the unitig into protein (up to six sequences),
    ## second: decompose into k-mers, save k-mers
    ## third: look for overlaps with query_kmers

    # cDBG records matching to queries:
    matching_cdbg = set()
    # total number of k-mers in cDBG:
    cdbg_kmers = 0

    db = sqlite3.connect(args.unitigs_db)
    for n, record in enumerate(search_utils.contigs_iter_sqlite(db)):
        if n and n % 1000 == 0:
            print(
                f"... searched {n} unitigs, found {len(matching_cdbg)} matching so far",
                end="\r",
                file=sys.stderr,
            )

        # translate into protein sequences
        seq = record.sequence
        record_hashes = translate_mh.seq_to_hashes(record.sequence)
        record_hashes = set(record_hashes)
        cdbg_kmers += len(record_hashes)

        # do we have an overlap with query??
        if record_hashes & query_kmers:
            # if so, save etc.
            matching_cdbg.add(record.name)

            # track sequence in DNA space...
            dna_mh.add_sequence(record.sequence)
            # ...and in protein space.
            prot_mh.add_sequence(record.sequence)
    n += 1

    print(f"loaded {n} query sequences (unitigs)", file=sys.stderr)
    print(f"total matches: {len(matching_cdbg)}", file=sys.stderr)

    dna_sig = sourmash.SourmashSignature(dna_mh, name="DNA")
    prot_sig = sourmash.SourmashSignature(prot_mh, name="prot")

    outsig = out_prefix + ".sig"
    with open(outsig, "wt") as fp:
        sourmash.save_signatures([dna_sig, prot_sig], fp)
    print(f"saved prot & dna signatures to '{outsig}'")

    outnodes = out_prefix + ".nodes.gz"
    with gzip.open(outnodes, "wt") as fp:
        id_list = "\n".join([str(x) for x in sorted(matching_cdbg)])
        print(id_list, file=fp)
    print(f"saved {len(matching_cdbg)} matching nodes to '{outnodes}'")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
