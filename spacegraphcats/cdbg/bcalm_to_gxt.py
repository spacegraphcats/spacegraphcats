#! /usr/bin/env python
"""
Convert a bcalm unitigs.fa output (a cDBG) into spacegraphcats files.

Outputs a GXT file (containing an undirected graph), a BGZF file
containing the sequences, and a .info.csv file containing
the BGZF offset, mean abundance, and length of each contig.

Also outputs sourmash k=31 scaled=1000 signatures for both input and
output files.
"""
import screed
import sys
import collections
import argparse
from spacegraphcats.utils.bgzf import bgzf
import logging
from typing import List, Dict, Set
import sourmash


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


def is_directed_path(x, sequences, neighbors, k):
    assert len(neighbors[x]) == 2,\
            ValueError("is_directed_path requires a degree 2 vertex")
    u, v = neighbors[x]
    seq_u = sequences[u]
    seq_v = sequences[v]
    if end_match(seq_u, seq_v, k, 'pp') or\
            end_match(seq_u, seq_v, k, 'ss'):
        return False
    else:
        return True


def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(comp[x] for x in reversed(seq))


def contract_neighbor(x, u, neighbors, sequences, mean_abunds, sizes, k):
    seq_x = sequences[x]
    seq_u = sequences[u]
    r = reverse_complement(seq_u)
    logging.debug("neighbor {}".format(u))
    logging.debug("seq_u: {}, {}".format(seq_u[:k-1], seq_u[1-k:]))
    logging.debug("   rc: {}, {}".format(r[:k-1], r[1-k:]))
    # check which direction the match occurs
    match, rc = end_match(seq_x, seq_u, k, 'sp')
    if match:
        if rc:
            sequences[x] += reverse_complement(seq_u)[k-1:]
        else:
            sequences[x] += seq_u[k-1:]
    else:
        # need to check whether reverse complement was used
        if seq_x[:k-1] == seq_u[1-k:]:
            sequences[x] = seq_u[:1-k] + seq_x
        else:
            r = reverse_complement(seq_u)
            assert seq_x[:k-1] == r[1-k:]
            sequences[x] = r[:1-k] + seq_x

    total_abunds = mean_abunds[x] * sizes[x] + \
        mean_abunds[u] * sizes[u]
    sizes[x] += sizes[u]
    mean_abunds[x] = total_abunds/sizes[x]
    # remove v from the graph by making u's other neighbor the
    # neighbor of x
    # there is no neighbor if u has degree 1
    if len(neighbors[u]) > 1:
        for y in neighbors[u]:
            if y != x:
                break
        else:
            msg = "{} doesn't have {} as a neighbor".format(u, x)
            raise ValueError(msg)
        neighbors[y].remove(u)
        neighbors[y].add(x)
        neighbors[x].add(y)
    else:
        y = None  # for debug purposes
    neighbors[x].remove(u)
    neighbors[u] = set()
    logging.debug(
        "removed {}, replacing it with {}, {}".format(u, x, y))


def contract_degree_two(non_pendants, neighbors, sequences, mean_abunds, sizes,
                        k):
    deg_2 = list()
    for v, N in sorted(neighbors.items()):   # do we need sorted here!?
        if v in non_pendants or len(N) == 0:
            continue
        u = list(N)[0]
        neighbors[u].remove(v)
        N.remove(u)
        if len(neighbors[u]) == 2:
            deg_2.append(u)
    logging.debug("len(deg_2): {}".format(len(deg_2)))
    for x in deg_2:
        if len(neighbors[x]) != 2:
            continue
        u, v = list(neighbors[x])
        seq_x = sequences[x]
        seq_u = sequences[u]
        seq_v = sequences[v]
        # if uxv doesn't form a directed path, we can't do anything
        if end_match(seq_u, seq_v, k, 'pp')[0] or\
                end_match(seq_u, seq_v, k, 'ss')[0]:
            continue
        logging.debug("analyzing {}".format(x))
        logging.debug("seq_x: {}, {}".format(seq_x[:k-1],
                      seq_x[1-k:]))
        # can only delete u or v if they have low degree and their
        # neighbors have a directed path
        u_deg = len(neighbors[u])
        if u_deg == 1 or (u_deg == 2 and
                          is_directed_path(u, sequences, neighbors, k)):
            contract_neighbor(x, u, neighbors, sequences, mean_abunds,
                              sizes, k)
            non_pendants.remove(u)
        v_deg = len(neighbors[v])
        if v_deg == 1 or (v_deg == 2 and
                          is_directed_path(v, sequences, neighbors, k)):
            contract_neighbor(x, v, neighbors, sequences, mean_abunds,
                              sizes, k)
            non_pendants.remove(v)


def read_bcalm(unitigs, debug, k):
    # track offsets, mean abunds, and # k-mers for each contig
    # track links between contig IDs
    neighbors = collections.defaultdict(set)
    mean_abunds = {}  # type: Dict[int, float]
    sizes = {}  # type: Dict[int, int]
    sequences = {}  # type: Dict[int, Text]

    # walk the input unitigs file, tracking links between contigs.
    print('reading unitigs from {}'.format(unitigs))
    for n, record in enumerate(screed.open(unitigs)):
        if n % 10000 == 0:
            print('...', n, file=sys.stderr, end='\r')

        name = record.name
        name_split = name.split()

        # note: contig_id may not be in order.
        contig_id = int(name_split[0])

        # track the various links
        links = [x for x in name_split[1:] if x.startswith('L:')]
        link_ids = [x.split(':')[2] for x in links]
        link_ids = [int(x) for x in link_ids if int(x) != contig_id]

        logging.debug('link_ids for {} are {}'.format(contig_id, link_ids))

        neighbors[contig_id].update(link_ids)

        # get mean abund
        abund = [x for x in name_split[1:] if x.startswith('km:')]
        assert len(abund) == 1, abund
        abund = abund[0].split(':')
        assert len(abund) == 3
        abund = float(abund[2])

        # record the abundance sequence and size
        mean_abunds[contig_id] = abund
        sequences[contig_id] = record.sequence
        sizes[contig_id] = len(record.sequence) - k + 1

    print('...read {} unitigs'.format(len(sequences)))

    fail = False
    for source in neighbors:
        for nbhd in neighbors[source]:
            if source not in neighbors[nbhd]:
                print('{} -> {}, but not {} -> {}'.format(source, nbhd,
                                                          nbhd, source))
                fail = True

    assert not fail

    return neighbors, sequences, mean_abunds, sizes


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
    parser.add_argument('--randomize', help='randomize cDBG order')
    args = parser.parse_args(argv)

    k = args.ksize

    trim = not args.pendants
    trim_cutoff = args.abundance
    unitigs = args.bcalm_unitigs
    debug = args.debug

    if args.debug:
        logging.basicConfig(filename='bcalm_to_gxt.log', filemode='w',
                            level=logging.DEBUG)
    else:
        logging.basicConfig(filename='bcalm_to_gxt.log', filemode='w',
                            level=logging.WARNING)

    logging.debug("starting bcalm_to_gxt run.")

    gxtfp = open(args.gxt_out, 'wt')
    contigsfp = bgzf.open(args.contigs_out, 'wb')
    info_filename = args.contigs_out + '.info.csv'
    info_fp = open(info_filename, 'wt')
    in_mh = sourmash.MinHash(0, 31, scaled=1000)
    out_mh = sourmash.MinHash(0, 31, scaled=1000)

    # load in the basic graph structure from the BCALM output file
    neighbors, sequences, mean_abunds, sizes = read_bcalm(unitigs, debug, k)

    # record input k-mers in a minhash
    for seq in sequences.values():
        in_mh.add_sequence(seq)

    # make order deterministic by reordering around min value of first, last,
    # and reverse complementing sequences appropriately
    print('reordering...')
    reordering = {}

    # first, put sequences in specific orientation
    sequence_list = []
    for key in neighbors:
        v = sequences[key]

        # pick lexicographically smaller of forward & reverse complement.
        v2 = screed.rc(v)
        if v > v2:
            v = v2
        sequence_list.append((v, key))
        del sequences[key]

    # sort all sequences:
    sequence_list.sort(reverse=True)
    if args.randomize:
        print('(!! randomizing order per --randomize !!)')
        random.shuffle(sequence_list)

    # ok, now remap all the things.
    remapping = {}
    new_sequences = {}

    # remap sequences
    new_key = 0
    while sequence_list:                  # consume while iterating
        sequence, old_key = sequence_list.pop()
        remapping[old_key] = new_key
        new_sequences[new_key] = sequence
        new_key += 1

    # remap other things
    new_neighbors = collections.defaultdict(set)
    for old_key, vv in neighbors.items():
        new_vv = [ remapping[v] for v in vv ]
        new_neighbors[remapping[old_key]] = set(new_vv)

    new_mean_abunds = {}
    for old_key, value in mean_abunds.items():
        new_mean_abunds[remapping[old_key]] = value

    new_sizes = {}
    for old_key, value in sizes.items():
        new_sizes[remapping[old_key]] = value

    assert len(sequences) == 0
    print('...done')

    sequences = new_sequences
    mean_abunds = new_mean_abunds
    sizes = new_sizes
    neighbors = new_neighbors

    # if we are removing pendants, we need to relabel the contigs so they are
    # consecutive integers starting from 0.  If not, we create dummy data
    # structures to make the interface the same elsewhere in the data
    if trim:
        print('removing pendants...')
        non_pendants = set(v for v, N in neighbors.items() if len(N) > 1 or
                           mean_abunds[v] > trim_cutoff)
        contract_degree_two(non_pendants, neighbors, sequences, mean_abunds,
                            sizes, k)
    else:
        non_pendants = list(neighbors.keys())
    aliases = {x: i for i, x in enumerate(sorted(non_pendants))}
    n = len(aliases)

    # write out sequences & compute offsets
    offsets = {}
    kv_list = sorted(aliases.items(), key=lambda x:x[1])
    for x, i in kv_list:
        offsets[x] = contigsfp.tell()
        contigsfp.write('>{}\n{}\n'.format(i, sequences[x]))
        out_mh.add_sequence(sequences[x])
    contigsfp.close()

    print('... done! {} unitigs'.format(n))

    # start the gxt file by writing the number of nodes (unitigs))
    gxtfp.write('{}\n'.format(n))

    # write out all of the links, in 'from to' format.
    n_edges = 0
    for v, N in sorted(neighbors.items()):
        for u in sorted(N):
            gxtfp.write('{} {}\n'.format(aliases[v], aliases[u]))
            n_edges += 1

    print('{} vertices, {} edges'.format(n, n_edges))

    info_fp.write('contig_id,offset,mean_abund,n_kmers\n')
    for v, i in aliases.items():
        info_fp.write('{},{},{:.3f},{}\n'.format(i, offsets[v],
                                                 mean_abunds[v],
                                                 sizes[v]))

    # output two sourmash signatures: one for input contigs, one for
    # output contigs.
    in_sig = sourmash.SourmashSignature(in_mh, filename=args.bcalm_unitigs)
    sourmash.save_signatures([ in_sig ],
                             open(args.bcalm_unitigs + '.sig', 'wt'))

    out_sig = sourmash.SourmashSignature(out_mh, filename=args.contigs_out)
    sourmash.save_signatures([ out_sig ],
                             open(args.contigs_out + '.sig', 'wt'))


if __name__ == '__main__':
    main(sys.argv[1:])
