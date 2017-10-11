#! /usr/bin/env python
import argparse
import khmer, khmer.utils
import screed


def main():
    p = argparse.ArgumentParser()
    p.add_argument('contigs')
    p.add_argument('reads')
    p.add_argument('reads_out')
    args = p.parse_args()

    print('consuming contigs:', args.contigs)
    kh = khmer.Nodegraph(21, 4e8, 4)
    kh.consume_seqfile(args.contigs)

    outfp = open(args.reads_out, 'wt')
    wrote = 0

    for n, record in enumerate(khmer.utils.clean_input_reads(screed.open(args.reads))):
        if n % 100000 == 0:
            print('...', n, wrote)

        if not kh.get_min_count(record.cleaned_seq):
            outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
            wrote += 1


if __name__ == '__main__':
    main()
