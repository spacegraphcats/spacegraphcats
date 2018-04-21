#! /usr/bin/env python
"""
Convert a bcalm unitigs.fa output (a cDBG) into spacegraphcats files.

Outputs a GXT file (containing an undirected graph), a BGZF file
containing the sequences, and a .info.csv file containing
the BGZF offset, mean abundance, and length of each contig.
"""
import screed
import sys
import collections
import argparse
from search.bgzf import bgzf

TRIM_CUTOFF = 1.1

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('bcalm_unitigs')
    parser.add_argument('gxt_out')
    parser.add_argument('contigs_out')
    parser.add_argument('-k', '--ksize', type=int, default=31)
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-P', '--pendants', action="store_true",
                        help="don't remove low abundance pendants")
    args = parser.parse_args(argv)

    ksize = args.ksize

    trim = not args.pendants

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
    sequences = {}

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

        if args.debug:
            print('link_ids for {} are {}'.format(contig_id, link_ids))

        link_d[contig_id].update(link_ids)

        # get mean abund
        abund = [ x for x in name_split[1:] if x.startswith('km:') ]
        assert len(abund) == 1, abund
        abund = abund[0].split(':')
        assert len(abund) == 3
        abund = float(abund[2])

        mean_abunds[contig_id] = abund

        # where are we in the output file?
        assert contig_id not in offsets

        sequences[contig_id] = record.sequence

        sizes[contig_id] = len(record.sequence) - ksize + 1

    # if we are removing pendants, we need to relabel the contigs so they are
    # consecutive integers starting from 0.  If not, we create dummy data
    # structures to make the interface the same elsewhere in the data
    if trim:
        non_pendants = [x for x, N in link_d.items() if len(N) > 1 and \
                        mean_abunds[x] > TRIM_CUTOFF]
    else:
        non_pendants = list(link_d.keys())
    aliases = {x: i for i, x in enumerate(non_pendants)}
    n = len(aliases)

    for x, i in aliases.items():
        contigsfp.write('>{}\n{}\n'.format(i, sequences[x]))
        offsets[x] = contigsfp.tell()
    contigsfp.close()

    print('... done! {} unitigs'.format(n))

    # start the gxt file by writing the number of nodes (unitigs))
    gxtfp.write('{}\n'.format(n))

    # write out all of the links, in 'from to' format.
    n_edges = 0
    for node, edgelist in link_d.items():
        if node not in aliases:
            continue
        for next_node in edgelist:
            if next_node not in aliases:
                continue
            gxtfp.write('{} {}\n'.format(aliases[node], aliases[next_node]))
            n_edges += 1

    print('{} vertices, {} edges'.format(n, n_edges))

    info_fp.write('contig_id,offset,mean_abund,n_kmers\n')
    for v, i in aliases.items():
        info_fp.write('{},{},{:.3f},{}\n'.format(i, offsets[v],
                                                 mean_abunds[v],
                                                 sizes[v]))


if __name__ == '__main__':
    main(sys.argv[1:])
