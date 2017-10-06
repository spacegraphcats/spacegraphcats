#! /usr/bin/env python
import sys
import khmer, khmer.utils
import screed

def kmers(seq, ksize):
    for i in range(len(seq) - ksize + 1):
        yield i, seq[i:i+ksize]



a = sys.argv[1]
b = sys.argv[2]

kh = khmer.Nodegraph(31, 1e8, 1)

for record in khmer.utils.clean_input_reads(screed.open(a)):
    kh.consume(record.cleaned_seq)

for record in khmer.utils.clean_input_reads(screed.open(b)):
    for pos, kmer in kmers(record.cleaned_seq, kh.ksize()):
        if kh.get(kmer) == 0:
            print(pos, kmer, kh.hash(kmer))
