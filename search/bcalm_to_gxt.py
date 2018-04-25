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

def find_thin_path(first, second, link_d):
    """
    Find the longest thin path whose first two vertices are (first, second.)
    """
    path = list()
    curr = second
    prev = first
    while len(link_d[curr]) == 2:
        path.append(curr)
        for next_v in link_d[curr]:
            if next_v != prev:
                break
        else:
            err = "{} has two neighbors, {}, but they are both"
                " {}".format(curr, link_d[curr], prev)
            raise ValueError(err)
        prev = curr
        curr = next_v
    return path


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('bcalm_unitigs')
    parser.add_argument('gxt_out')
    parser.add_argument('contigs_out')
    parser.add_argument('-k', '--ksize', type=int, default=31)
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-P', '--pendants', action="store_true",
                        help="don't remove low abundance pendants")
    parser.add_argument('-a', '--abundance', nargs='?', type=float,
                        default=1.1)
    args = parser.parse_args(argv)

    ksize = args.ksize

    trim = not args.pendants
    trim_cutoff = args.abundance

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


        # record the abundance sequence and size
        mean_abunds[contig_id] = abund
        sequences[contig_id] = record.sequence
        sizes[contig_id] = len(record.sequence) - ksize + 1

    # if we are removing pendants, we need to relabel the contigs so they are
    # consecutive integers starting from 0.  If not, we create dummy data
    # structures to make the interface the same elsewhere in the data
    if trim:
        non_pendants = set(v for v, N in link_d.items() if len(N) > 1 or \
                        mean_abunds[v] > trim_cutoff)
        for v, N in link_d.items():
            if v in non_pendants or len(N) == 0:
                continue
            u = N[0]
            link_d[u].remove(v)

        # pendant removal may have made some vertices with degree 2.  If this creates a long path whose vertices all have degree 2 (a thin path), we should replace the thin path by a single vertex
        thin_paths = set()  # vertices in a thin path
        for x in non_pendants:
            if x in thin_paths or len(link_d[x]) != 2:
                continue
            u, v = link_d[x]
            # start our searches with vertices that have a neighbor of high degree and a neighbor of degree two.  u and v are thus far indistinguishable so we have to try both possibilities (v has high degree vs. u has high degree)
            if len(link_d[u]) != 2:
                path = find_thin_path(x, v, link_d)
            elif len(link_d[v]) != 2:
                path = find_thin_path(x, u, link_d)
            else:
                path = list()
            thin_paths.update(path)
            # update sequences
            for w in path:
                # decide whether the suffix of w is a prefix of x or vice versa
                if sequences[w][:ksize-1] == sequences[x][1-ksize:]:
                    # prefix of w is a suffix of x
                    sequences[x] += sequences[w][ksize-1:]
                elif sequences[x][:ksize-1] == sequences[w][1-ksize:]:
                    # prefix of x is a suffix of w
                    sequences[x] = sequences[w][:1-ksize] + sequences[x]
                else:
                    # there's an error
                    err = "{} and {} don't have prefix-suffix overlap".format(
                        sequences[w], sequences[x])
                    raise ValueError(err)
                # update sizes and mean_abunds
                total_abunds = mean_abunds[x] * sizes[x] + \
                    mean_abunds[w] * (sizes[w] - ksize + 1)
                sizes[x] += sizes[w] - ksize + 1
                mean_abunds[x] = total_abunds/sizes[x]

            # update neighbors
            if len(path) > 0:
                for y in link_d[w]:
                    if len(link_d[y]) != 2:
                        break
                link_d[y].remove(w)
                link_d[x].remove(path[0])
                link_d[x].append(y)
                link_d[y].append(x)
        non_pendants -= thin_paths
    else:
        non_pendants = list(link_d.keys())
    aliases = {x: i for i, x in enumerate(non_pendants)}
    n = len(aliases)

    for x, i in aliases.items():
        offsets[x] = contigsfp.tell()
        contigsfp.write('>{}\n{}\n'.format(i, sequences[x]))
    contigsfp.close()

    print('... done! {} unitigs'.format(n))

    # start the gxt file by writing the number of nodes (unitigs))
    gxtfp.write('{}\n'.format(n))

    # write out all of the links, in 'from to' format.
    n_edges = 0
    for node, edgelist in link_d.items():
        # if node not in aliases:
        #     continue
        for next_node in edgelist:
            # if next_node not in aliases:
            #     continue
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
