#! /usr/bin/env python
import screed
import sys
import collections
import argparse
from search.bgzf import bgzf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bcalm_unitigs')
    parser.add_argument('gxt_out')
    parser.add_argument('contigs_out')
    parser.add_argument('-k', '--ksize', type=int, default=31)
    args = parser.parse_args()

    ksize = args.ksize

    # track links between contig IDs
    link_d = collections.defaultdict(set)

    gxtfp = open(args.gxt_out, 'wt')
    contigsfp = bgzf.open(args.contigs_out, 'wb')
    info_filename = args.contigs_out + '.info.csv'
    info_fp = open(info_filename, 'wt')

    # track offsets, mean abunds, and # k-mers for each contig
    offsets = {}
    mean_abunds = {}
    sizes = {}

    # walk the input unitigs file, tracking links between contigs and
    # writing them to contigs_out.
    max_contig_id = 0
    print('reading unitigs from {}'.format(args.bcalm_unitigs))
    for n, record in enumerate(screed.open(args.bcalm_unitigs)):
        if n % 10000 == 0:
            print('...', n, file=sys.stderr, end='\r')

        name = record.name
        name_split = name.split()

        # note: contig_id may not be in order.
        contig_id = int(name_split[0])

        # track the various links
        links = [ x for x in name_split[1:] if x.startswith('L:') ]
        link_ids = [ x.split(':')[2] for x in links ]
        link_ids = [ int(x) for x in link_ids ]

        link_d[contig_id].update(link_ids)

        # get mean abund
        abund = [ x for x in name_split[1:] if x.startswith('km:') ]
        assert len(abund) == 1, abund
        abund = float(abund[0].split(':')[2])

        mean_abunds[contig_id] = abund

        # where are we in the output file?
        assert contig_id not in offsets
        offsets[contig_id] = contigsfp.tell()

        # write out contig in boring format
        contigsfp.write('>{}\n{}\n'.format(contig_id, record.sequence))

        max_contig_id = max(contig_id, max_contig_id)

        sizes[contig_id] = len(record.sequence) - ksize + 1

    contigsfp.close()

    print('... done! {} unitigs'.format(n))
    assert max_contig_id == n
    assert max(link_d.keys()) <= n

    # start the gxt file by writing the number of nodes (unitigs))
    gxtfp.write('{}\n'.format(max(link_d.keys()) + 1))

    # write out all of the links, in 'from to' format.
    n_edges = 0
    for k, v in link_d.items():
        for vv in v:
            assert k <= n
            assert vv <= n
            gxtfp.write('{} {}\n'.format(k, vv))
            n_edges += 1

    print('{} vertices, {} edges'.format(n, n_edges))

    info_fp.write('contig_id,offset,mean_abund,n_kmers\n')
    for contig_id in range(max_contig_id + 1):
        info_fp.write('{},{},{:.3f},{}\n'.format(contig_id, offsets[contig_id],
                                                 mean_abunds[contig_id],
                                                 sizes[contig_id]))


if __name__ == '__main__':
    main()
