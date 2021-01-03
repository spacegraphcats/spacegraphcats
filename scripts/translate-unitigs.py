#! /usr/bin/env python
import sys
import screed
import argparse

dna_to_aa={'TTT':'F','TTC':'F', 'TTA':'L','TTG':'L',
                'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                'TAT':'Y','TAC':'Y', 'TAA':'*','TAG':'*','TGA':'*',
                'TGT':'C','TGC':'C', 'TGG':'W',
                'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
                'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                'CAT':'H','CAC':'H', 'CAA':'Q','CAG':'Q',
                'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                'ATT':'I','ATC':'I','ATA':'I', 'ATG':'M',
                'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                'AAT':'N','AAC':'N', 'AAA':'K','AAG':'K',
                'AGT':'S','AGC':'S', 'AGA':'R','AGG':'R',
                'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
                'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                'GAT':'D','GAC':'D', 'GAA':'E','GAG':'E',
                'GGT':'G','GGC':'G','GGA':'G','GGG':'G'}


__complementTranslation = { "A": "T", "C": "G", "G": "C", "T": "A", "N": "N" }
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
        yield dna_to_aa.get(seq[i:i+3], "X")


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
    p.add_argument('unitigs')
    p.add_argument('-k', '--ksize', type=int, default=10)
    args = p.parse_args()

    dna_ksize = args.ksize * 3

    for n, record in enumerate(screed.open(args.unitigs)):
        all_kmers = set()
        if n and n % 1000 == 0:
            print(f'... {n}', end='\r', file=sys.stderr)

        seq = record.sequence
        seqlen = len(seq)
        if seqlen < dna_ksize:
            continue
        elif seqlen == dna_ksize:
            prots = list(translate(seq, 1))
        elif seqlen == dna_ksize + 1:
            prots = list(translate(seq, 2))
        else:
            prots = list(translate(seq, 3))

        for prot in prots:
            these_kmers = set(kmers(prot, args.ksize))
            all_kmers.update(these_kmers)

        print(seqlen, len(all_kmers))

    return 0


if __name__ == '__main__':
    sys.exit(main())
