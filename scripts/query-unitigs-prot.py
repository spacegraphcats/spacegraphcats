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


dna_to_aa = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
             'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
             'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGA': '*',
             'TGT': 'C', 'TGC': 'C', 'TGG': 'W',
             'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
             'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
             'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
             'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
             'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
             'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
             'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
             'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
             'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
             'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
             'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
             'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

__complementTranslation = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def complement(s):
    """
    Return complement of 's'.
    """
    c = "".join(__complementTranslation[n] for n in s)
    return c


def reverse(s):
    """
    Return reverse of 's'.
    """
    r = "".join(reversed(s))
    return r


def peptides(seq, start):
    for i in range(start, len(seq), 3):
        yield dna_to_aa.get(seq[i:i + 3], "X")


def translate(seq, n_frames):
    for i in range(n_frames):
        pep = peptides(seq, i)
        yield "".join(pep)

    revcomp = reverse(complement((seq)))
    for i in range(n_frames):
        pep = peptides(revcomp, i)
        yield "".join(pep)


def kmers(seq, k):
    for start in range(len(seq) - k + 1):
        yield seq[start:start + k]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query')
    p.add_argument('unitigs_db')
    p.add_argument('-k', '--ksize', type=int, default=10,
                   help='protein ksize')
    p.add_argument('--out-prefix')
    args = p.parse_args()

    out_prefix = args.out_prefix
    if not out_prefix:
        out_prefix = os.path.basename(args.query)

    prot_ksize = args.ksize
    dna_ksize = args.ksize * 3

    dna_mh = sourmash.MinHash(n=0, scaled=1000, ksize=dna_ksize)
    prot_mh = sourmash.MinHash(n=0, scaled=100, ksize=dna_ksize,
                               is_protein=True)

    ## query: for all of the query sequences (which are protein),
    ## translate them into k-mers at protein ksize.
    query_kmers = set()
    for n, record in enumerate(screed.open(args.query)):
        if n and n % 1000 == 0:
            print(f'... {n} {len(query_kmers)} (query)', end='\r', file=sys.stderr)
        these_kmers = kmers(record.sequence, prot_ksize)
        query_kmers.update(these_kmers)

    print(f'loaded {n} query sequences {len(query_kmers)} prot kmers (query)', file=sys.stderr)

    ## now unitigs... for all of the unitigs (which are DNA),
    ## first: translate the unitig into protein (up to six sequences),
    ## second: decompose into k-mers, save k-mers
    ## third: look for overlaps with query_kmers

    matching_cdbg = set()
    cdbg_kmers = 0

    db = sqlite3.connect(args.unitigs_db)
    for n, record in enumerate(search_utils.contigs_iter_sqlite(db)):
        if n and n % 1000 == 0:
            print(f'... searched {n} unitigs, found {len(matching_cdbg)} matching so far', end='\r', file=sys.stderr)

        # translate into protein sequences
        seq = record.sequence
        seqlen = len(seq)
        if seqlen < dna_ksize:
            continue
        elif seqlen == dna_ksize:
            prots = list(translate(seq, 1))
            assert len(prots) == 2
        elif seqlen == dna_ksize + 1:
            prots = list(translate(seq, 2))
            assert len(prots) == 4
        else:
            prots = list(translate(seq, 3))
            assert len(prots) == 6

        # convert into k-mers
        record_kmers = set()
        for prot in prots:
            these_kmers = set(kmers(prot, prot_ksize))
            record_kmers.update(these_kmers)

        cdbg_kmers += len(record_kmers)

        # do we have an overlap with query??
        if record_kmers & query_kmers:
            # if so, save etc.
            matching_cdbg.add(record.name)
            dna_mh.add_sequence(record.sequence)
            for prot in prots:
                prot_mh.add_protein(prot)

    print(f'loaded {n} query sequences (unitigs)', file=sys.stderr)
    print(f'total matches: {len(matching_cdbg)}', file=sys.stderr)

    dna_sig = sourmash.SourmashSignature(dna_mh, name="DNA")
    prot_sig = sourmash.SourmashSignature(prot_mh, name="prot")

    outsig = out_prefix + '.sig'
    with open(outsig, 'wt') as fp:
        sourmash.save_signatures([dna_sig, prot_sig], fp)
    print(f"saved prot & dna signatures to '{outsig}'")

    outnodes = out_prefix + '.nodes.gz'
    with gzip.open(outnodes, "wt") as fp:
        id_list = "\n".join([str(x) for x in sorted(matching_cdbg)])
        print(id_list, file=fp)
    print(f"saved {len(matching_cdbg)} matching nodes to '{outnodes}'")

    return 0


if __name__ == '__main__':
    sys.exit(main())
