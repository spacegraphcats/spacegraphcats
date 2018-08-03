#! /usr/bin/env python
"""
Look for catlas nodes that have few k-mers but many cDBG nodes underneath
them, as a likely sign of strain variation.

EXPERIMENTAL.
"""
import argparse
import os
import sys
import math
import sourmash_lib

import screed
from .catlas import CAtlas


def main(args=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('output')
    p.add_argument('--minsize', type=float, default=100)
    p.add_argument('--maxsize', type=float, default=10000)
    p.add_argument('--keep-fraction', type=float, default=0.1)
    p.add_argument('-k', '--ksize', default=31, type=int,
                   help='k-mer size (default: 31)')
    args = p.parse_args(args)

    print('minsize: {:g}'.format(args.minsize))
    print('maxsize: {:g}'.format(args.maxsize))

    catlas_file = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')
    sizefile = os.path.join(args.catlas_prefix, 'contigs.fa.gz.info.csv')

    # load catlas DAG
    catlas = CAtlas(catlas_file, domfile=domfile, sizefile=sizefile)
    print('loaded {} nodes from catlas {}'.format(len(catlas), catlas_file))
    print('loaded {} layer 1 catlas nodes'.format(len(catlas.layer1_to_cdbg)))

    # calculate the cDBG shadow sizes for each catlas node.
    print('decorating catlas with shadow size info.')
    catlas.decorate_with_shadow_sizes()

    # ok, the real work: look at articulation of cDBG graph.

    # find highest nodes with kmer size less than given max_size
    def find_terminal_nodes(node_id, max_size):
        node_list = set()
        for sub_id in catlas.children[node_id]:
            # shadow size
            size = catlas.kmer_sizes[sub_id]

            if size < max_size:
                node_list.add(sub_id)
            else:
                children = find_terminal_nodes(sub_id, max_size)
                node_list.update(children)

        return node_list

    print('finding terminal nodes for {}.'.format(args.maxsize))

    terminal = find_terminal_nodes(catlas.root, args.maxsize)
    print('...got {}'.format(len(terminal)))
    terminal = {n for n in terminal if catlas.kmer_sizes[n] > args.minsize}
    print('...down to {} between {} and {} in size.'.format(len(terminal),
                                                            args.minsize,
                                                            args.maxsize))

    # now, go through and calculate ratios
    x = []
    for node_id in terminal:
        # calculate: how many k-mers per cDBG node?
        kmer_size = catlas.kmer_sizes[node_id]
        shadow_size = catlas.shadow_sizes[node_id]

        ratio = math.log(kmer_size, 2) - math.log(shadow_size, 2)

        # track basic info
        x.append((ratio, node_id, shadow_size, kmer_size))

    print('terminal node stats for maxsize: {:g}'.format(args.maxsize))
    print('n tnodes:', len(terminal))
    print('total k-mers:', catlas.kmer_sizes[catlas.root])

    x.sort(reverse=True)
    for (k, v, a, b) in x[:10]:
        print('ratio: {:.3f}'.format(2**k), '/ shadow size:', a, '/ kmers:', b)
    print('... eliding {} nodes'.format(len(x) - 20))
    for (k, v, a, b) in x[-10:]:
        print('ratio: {:.3f}'.format(2**k), '/ shadow size:', a, '/ kmers:', b)

    # keep the last keep-fraction (default 10%) for examination
    keep_sum_kmer = args.keep_fraction * catlas.kmer_sizes[catlas.root]
    sofar = 0
    keep_terminal = set()
    for (k, v, a, b) in reversed(x):
        sofar += b
        if sofar > keep_sum_kmer:
            break
        keep_terminal.add(v)

    print('keeping last {} k-mers worth of nodes for'
          'examination.'.format(sofar))

    # build cDBG shadow ID list.
    cdbg_shadow = catlas.shadow(keep_terminal)

    # extract contigs
    print('extracting contigs & building a sourmash signature')
    contigs = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

    # track results as signature
    contigs_mh = sourmash_lib.MinHash(n=0, ksize=args.ksize, scaled=1000)

    total_bp = 0
    total_seqs = 0

    outfp = open(args.output, 'wt')
    for n, record in enumerate(screed.open(contigs)):
        if n and n % 10000 == 0:
            offset_f = total_seqs / len(cdbg_shadow)
            print('...at n {} ({:.1f}% of shadow)'.format(total_seqs,
                  offset_f * 100),
                  end='\r')

        # contig names == cDBG IDs
        contig_id = int(record.name)
        if contig_id not in cdbg_shadow:
            continue

        outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
        contigs_mh.add_sequence(record.sequence)

        # track retrieved sequences in a minhash
        total_bp += len(record.sequence)
        total_seqs += 1

    # done - got all contigs!
    print('')
    print('fetched {} contigs, {} bp.'.format(total_seqs, total_bp))

    print('wrote contigs to {}'.format(args.output))
    with open(args.output + '.sig', 'wt') as fp:
        ss = sourmash_lib.SourmashSignature(contigs_mh)
        sourmash_lib.save_signatures([ss], fp)


if __name__ == '__main__':
    main()
