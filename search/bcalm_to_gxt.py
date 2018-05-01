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
import logging

def find_thin_path(first, second, link_d):
    """
    Find the longest thin path whose first two vertices are (first, second.)
    """
    path = list()
    curr = second
    prev = first
    while len(link_d[curr]) == 2:
        path.append(curr)
        N = link_d[curr]
        logging.debug("{} has degree 2 with neighbors: {}".format(curr, N))
        for next_v in N:
            if next_v != prev:
                break
        else:
            err = "{} has two neighbors, {}, but they are both {}".format(
                curr, N, prev)
            raise ValueError(err)
        prev = curr
        curr = next_v
    logging.debug("{} has degree {}".format(curr, len(link_d[curr]))) 
    return path


def end_match(s, t, k, direction='sp'):
    """
    Compares the first or last k-1 bases of strings s,t for a match.
    The direction argument is a two character string whose characters are 'p'
    or 's', where the first character indicates whether to use the prefix or
    suffix of s for comparison and the second character indicates likewise for
    t.
    Returns a boolean 2-tuple whose first entry indicates whether there is a
    match and whose second entry indicates whether the reverse complement was
    necessary for a match.  If there is no match, the second entry will always
    be False.
    """
    r = reverse_complement(t)
    if direction not in ['pp', 'ss', 'ps', 'sp']:
        raise ValueError("Valid directions are 'pp', 'ss', 'ps', 'sp'")
    if direction[0] == 'p':
        s_end = s[:k-1]
    else:
        s_end = s[1-k:]
    if direction[1] == 'p':
        t_end = t[:k-1]
        r_end = r[:k-1]
    else:
        t_end = t[1-k:]
        r_end = r[1-k:]
    if s_end == t_end:
        return (True, False)
    elif s_end == r_end:
        return (True, True)
    else:
        return (False, False)


def is_directed_path(x, sequences, neighbors, ksize):
    assert len(neighbors[x]) == 2,\
            ValueError("is_directed_path requires a degree 2 vertex")
    u, v = neighbors[x]
    seq_u = sequences[u]
    seq_v = sequences[v]
    if end_match(seq_u, seq_v, ksize, 'pp') or\
            end_match(seq_u, seq_v, ksize, 'ss'):
        return False
    else:
        return True


def reverse_complement(seq):
    comp = {'A': 'T', 'T':'A', 'C': 'G', 'G': 'C'}
    return "".join(comp[x] for x in reversed(seq))


def main(argv):
    logging.basicConfig(filename='bcalm_to_gxt.log', filemode='w',
                        level=logging.DEBUG)
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
        link_ids = [ int(x) for x in link_ids if int(x) != contig_id]

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
        dead_neighbors = collections.defaultdict(list)
        non_pendants = set(v for v, N in link_d.items() if len(N) > 1 or \
                        mean_abunds[v] > trim_cutoff)
        deg_2 = list()
        for v, N in link_d.items():
            if v in non_pendants or len(N) == 0:
                continue
            u = list(N)[0]
            link_d[u].remove(v)
            N.remove(u)
            pre = sequences[v][:ksize-1]
            suf = sequences[v][1-ksize:]
            dead_neighbors[u].append((v, pre, suf))
            if len(link_d[u]) == 2:
                deg_2.append(u)
        logging.debug("len(deg_2): {}".format(len(deg_2)))
        for x in deg_2:
            if len(link_d[x]) != 2:
                continue
            u, v = list(link_d[x])
            seq_x = sequences[x]
            seq_u = sequences[u]
            seq_v = sequences[v]
            # if uxv doesn't form a directed path, we can't do anything
            if end_match(seq_u, seq_v, ksize, 'pp')[0] or\
                    end_match(seq_u, seq_v, ksize, 'ss')[0]:
                continue
            logging.debug("analyzing {}".format(x))
            logging.debug("seq_x: {}, {}".format(seq_x[:ksize-1], seq_x[1-ksize:]))
            # can only delete u or v if they have low degree and their
            # neighbors have a directed path
            u_deg = len(link_d[u])
            v_deg = len(link_d[v])
            if u_deg == 1 or (u_deg == 2 and\
                    is_directed_path(u, sequences, link_d, ksize)):
                r = reverse_complement(seq_u)
                logging.debug("neighbor {}".format(u))
                logging.debug("seq_u: {}, {}".format(seq_u[:ksize-1], seq_v[1-ksize:]))
                logging.debug("   rc: {}, {}".format(r[:ksize-1], r[1-ksize:]))
                # check which direction the match occurs
                match, rc = end_match(seq_x, seq_u, ksize, 'sp')
                if match:
                    if rc:
                        sequences[x] += reverse_complement(seq_u)[ksize-1:]
                    else:
                        sequences[x] += seq_u[ksize-1:]
                else:
                    # need to check whether reverse complement was used
                    if seq_x[:ksize-1] == seq_u[1-ksize:]:
                        sequences[x] = seq_u[:1-ksize] + seq_x
                    else:
                        r = reverse_complement(seq_u)
                        assert seq_x[:ksize-1] == r[1-ksize:]
                        sequences[x] = r[:1-ksize] + seq_x

                total_abunds = mean_abunds[x] * sizes[x] + \
                    mean_abunds[u] * sizes[u] 
                sizes[x] += sizes[u]
                mean_abunds[x] = total_abunds/sizes[x]
                # remove v from the graph by making u's other neighbor the
                # neighbor of x
                # there is no neighbor if u has degree 1
                if len(link_d[u]) > 1:
                    for y in link_d[u]:
                        if y != x:
                            break
                    else:
                        msg = "{} doesn't have {} as a neighbor".format(v, x)
                        raise ValueError(msg)
                    link_d[y].remove(u)
                    link_d[y].add(x)
                    link_d[x].add(y)
                else:
                    y = None  # for debug purposes
                link_d[x].remove(u)
                link_d[u] = set()
                non_pendants.remove(u)
                logging.debug("removed {}, replacing it with {}, {}".format(u, x, y))

            seq_x = sequences[x]
            if v_deg == 1 or (v_deg == 2 and\
                    is_directed_path(v, sequences, link_d, ksize)):
                r = reverse_complement(seq_v)
                logging.debug("neighbor {}".format(v))
                logging.debug("seq_v: {}, {}".format(seq_v[:ksize-1], seq_v[1-ksize:]))
                logging.debug("   rc: {}, {}".format(r[:ksize-1], r[1-ksize:]))
                # check which direction the match occurs
                match, rc = end_match(seq_x, seq_v, ksize, 'sp')
                logging.debug("sp matching: {}, {}".format(match, rc))
                if match:
                    if rc:
                        sequences[x] += reverse_complement(seq_v)[ksize-1:]
                    else:
                        sequences[x] += seq_v[ksize-1:]
                else:
                    logging.debug("ps matching: {}".format(end_match(seq_x, seq_v, ksize, 'ps')))
                    # need to check whether reverse complement was used
                    if seq_x[:ksize-1] == seq_v[1-ksize:]:
                        sequences[x] = seq_v[:1-ksize] + seq_x
                    else:
                        r = reverse_complement(seq_v)
                        assert seq_x[:ksize-1] == r[1-ksize:]
                        sequences[x] = r[:1-ksize] + seq_x

                total_abunds = mean_abunds[x] * sizes[x] + \
                    mean_abunds[v] * sizes[v]
                sizes[x] += sizes[v]
                mean_abunds[x] = total_abunds/sizes[x]
                # remove v from the graph by making u's other neighbor the
                # neighbor of x
                # there is no neighbor if v has degree 1
                if len(link_d[v]) > 1:
                    for y in link_d[v]:
                        if y != x:
                            break
                    else:
                        msg = "{} doesn't have {} as a neighbor".format(v, x)
                        raise ValueError(msg)
                    link_d[y].remove(v)
                    link_d[y].add(x)
                    link_d[x].add(y)
                else:
                    y = None  # for debug purposes
                link_d[x].remove(v)
                link_d[v] = set()       
                non_pendants.remove(v)
                logging.debug("removing {}, replacing with {}, {}".format(v, x, y))
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
    for v, edgelist in link_d.items():
        # if node not in aliases:
        #     continue
        for u in edgelist:
            try:
                u_id = aliases[u]
            except KeyError as e:
                print("{} {}".format(v, edgelist))
                raise e
            gxtfp.write('{} {}\n'.format(aliases[v], aliases[u]))
            n_edges += 1

    print('{} vertices, {} edges'.format(n, n_edges))

    info_fp.write('contig_id,offset,mean_abund,n_kmers\n')
    for v, i in aliases.items():
        info_fp.write('{},{},{:.3f},{}\n'.format(i, offsets[v],
                                                 mean_abunds[v],
                                                 sizes[v]))


if __name__ == '__main__':
    main(sys.argv[1:])
