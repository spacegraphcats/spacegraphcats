#! /usr/bin/env python3
from __future__ import print_function

import sys
import khmer
import screed
import argparse
from collections import OrderedDict
import os, os.path
from spacegraphcats import graph_parser


# graph settings
DEFAULT_KSIZE=31
NODEGRAPH_SIZE=8e8

# minhash settings
MH_SIZE_DIVISOR=50
MH_MIN_SIZE=5

class Pathfinder(object):
    "Track segment IDs, adjacency lists, and MinHashes"
    def __init__(self, ksize, node_offset=0):
        self.ksize = ksize

        self.node_counter = 1 + node_offset
        self.nodes = {}                       # node IDs (int) to size
        self.nodes_to_kmers = {}              # node IDs (int) to kmers
        self.kmers_to_nodes = {}              # kmers to node IDs
        self.adjacencies = {}                 # node to node
        self.hashdict = OrderedDict()         # node to MinHash ID list
        self.labels = {}                      # nodes to set of labels

    def new_hdn(self, kmer):
        if kmer in self.kmers_to_nodes:
            return self.kmers_to_nodes[kmer]

        this_id = self.node_counter
        self.node_counter += 1

        self.nodes[this_id] = 1
        self.nodes_to_kmers[this_id] = kmer
        self.kmers_to_nodes[kmer] = this_id

        return this_id

    def new_linear_node(self, visited, size):
        node_id = self.node_counter
        self.node_counter += 1
        self.nodes[node_id] = size

        kmer = min(visited)               # represent linear nodes by min(hash)
        self.nodes_to_kmers[node_id] = kmer
        self.kmers_to_nodes[kmer] = node_id

        return node_id

    def add_adjacency(self, node_id, adj):
        node_id, adj = min(node_id, adj), max(node_id, adj)
        
        x = self.adjacencies.get(node_id, set())
        x.add(adj)
        self.adjacencies[node_id] = x

    def add_label(self, kmer, label):
        x = self.labels.get(kmer, set())
        x.add(label)
        self.labels[kmer] = x


def traverse_and_mark_linear_paths(graph, nk, stop_bf, pathy, degree_nodes):
    size, adj_kmers, visited = graph.traverse_linear_path(nk, degree_nodes,
                                                          stop_bf)
    if not size:
        return

    # get a ID
    path_id = pathy.new_linear_node(visited, size)

    # for all adjacencies, add.
    for kmer in adj_kmers:
        adj_node_id = pathy.kmers_to_nodes.get(kmer)
        pathy.add_adjacency(path_id, adj_node_id)

    # next, calculate minhash from visited k-mers
    v = [ khmer.reverse_hash(i, graph.ksize()) for i in visited ]
    mh_size = max(len(visited) // MH_SIZE_DIVISOR, MH_MIN_SIZE)

    mh = khmer.MinHash(mh_size, graph.ksize())
    for kmer in v:
        mh.add_sequence(kmer)

    # save minhash info
    pathy.hashdict[path_id] = mh.get_mins()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('seqfiles', nargs='+')
    p.add_argument('-o', '--output', default=None)
    p.add_argument('-k', '--ksize', default=DEFAULT_KSIZE, type=int)
    p.add_argument('-x', '--tablesize', default=NODEGRAPH_SIZE,
                            type=float)
    p.add_argument('--force', action='store_true')
    p.add_argument('--label', action='store_true')
    p.add_argument('--label-offset', type=int, default=0, help='for debug')
    p.add_argument('--node-offset', type=int, default=0, help='for debug')
    p.add_argument('--label-linear-segments', action='store_true')
    p.add_argument('--no-label-hdn', action='store_true')
    args = p.parse_args()

    assert args.ksize % 2, "ksize must be odd"
    if args.label_linear_segments or args.no_label_hdn:
        assert args.label

    output_dir = args.output
    if not output_dir:
        if len(args.seqfiles) > 1:
            print('** please specify an output directory with -o',
                  file=sys.stderr)
            sys.exit(-1)

        output_dir = os.path.basename(args.seqfiles[0])
        if output_dir.endswith('.fa'):
            output_dir = output_dir[:-3]

    gxtfile = os.path.basename(output_dir) + '.gxt'
    gxtfile = os.path.join(output_dir, gxtfile)
    mxtfile = os.path.basename(output_dir) + '.mxt'
    mxtfile = os.path.join(output_dir, mxtfile)

    print('')
    print('placing output in directory:', output_dir)
    print('gxt will be:', gxtfile)
    print('mxt will be:', mxtfile)
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        print('(note: directory already exists)')

    print('')
    print('building graphs and loading files')

    # Create graph, and two stop bloom filters - one for loading, one for
    # traversing. Create them all here so that we can error out quickly
    # if memory is a problem.

    graph = khmer.Nodegraph(args.ksize, args.tablesize, 2)
    stop_bf = khmer.Nodegraph(args.ksize, args.tablesize, 2)
    n = 0

    # load in all of the input sequences, one file at a time.
    for seqfile in args.seqfiles:
        for record in screed.open(seqfile):
            if len(record.sequence) < args.ksize: continue
            n += 1
            if n % 10000 == 0:
                print('...', seqfile, n)
            graph.consume(record.sequence)

    # complain if too small set of graphs was used.
    fp_rate = khmer.calc_expected_collisions(graph,
                                             args.force, max_false_pos=.05)

    # initialize the object that will track information for us.
    pathy = Pathfinder(args.ksize, args.node_offset)

    print('finding high degree nodes')
    if args.label and not args.no_label_hdn:
        print('(and labeling them, per request)')
    degree_nodes = khmer.HashSet(args.ksize)
    n = args.label_offset
    for seqfile in args.seqfiles:
        for record in screed.open(seqfile):
            if len(record.sequence) < args.ksize: continue
            n += 1
            if n % 10000 == 0:
                print('...2', seqfile, n)
            # walk across sequences, find all high degree nodes,
            # name them and cherish them.
            these_hdn = graph.find_high_degree_nodes(record.sequence)
            degree_nodes += these_hdn
            if args.label and not args.no_label_hdn:
                for kmer in these_hdn:
                    pathy.add_label(kmer, n)

    # get all of the degree > 2 kmers and give them IDs.
    for kmer in degree_nodes:
        pathy.new_hdn(kmer)

    print('traversing linear segments from', len(degree_nodes), 'nodes')

    # now traverse from each high degree nodes into all neighboring nodes,
    # seeking adjacencies.  if neighbor is high degree node, add it to
    # adjacencies; if neighbor is not, then traverse the linear path.  also
    # track minhashes while we're at it.
    for n, k in enumerate(degree_nodes):
        if n % 10000 == 0:
            print('...', n, 'of', len(degree_nodes))

        # retrieve the node ID of the primary segment.
        k_id = pathy.kmers_to_nodes[k]

        # add its minhash value.
        k_str = khmer.reverse_hash(k, graph.ksize())
        mh = khmer.MinHash(1, graph.ksize())
        mh.add_sequence(k_str)
        pathy.hashdict[k_id] = mh.get_mins()

        # find all the neighbors of this high-degree node.
        nbh = graph.neighbors(k)
        for nk in nbh:
            # neighbor is high degree? fine, mark its adjacencies.
            if nk in degree_nodes:
                nk_id = pathy.kmers_to_nodes[nk]
                pathy.add_adjacency(k_id, nk_id)
            else:
                # linear! walk it.
                traverse_and_mark_linear_paths(graph, nk, stop_bf, pathy,
                                               degree_nodes)

    print(len(pathy.nodes), 'segments, containing',
              sum(pathy.nodes.values()), 'nodes')

    if args.label and args.label_linear_segments:
        print('...doing labeling of linear segments by request.')
        n = args.label_offset
        for seqfile in args.seqfiles:
            for record in screed.open(seqfile):
                if len(record.sequence) < args.ksize: continue

                all_kmers = graph.get_kmer_hashes_as_hashset(record.sequence)
                n += 1

                for kmer, path_id in pathy.kmers_to_nodes.items():
                    if kmer in all_kmers:
                        pathy.add_label(kmer, n)

    # save to GXT/MXT.
    print('saving gxtfile', gxtfile)

    all_labels = set()
    label_counts = {}
    with open(gxtfile, 'w') as fp:
        w = graph_parser.Writer(fp, ['labels'], [])

        for k, v in pathy.nodes.items():
            kmer = pathy.nodes_to_kmers.get(k)
            l = ""
            if kmer:
                labels = pathy.labels.get(kmer, "")
                if labels:
                    for x in labels:
                        label_counts[x] = label_counts.get(x, 0) + 1
                    all_labels.update(labels)
                    l = " ".join(map(str, labels))
            w.add_vertex(k, v, [l])

        for k, v in pathy.adjacencies.items():
            for edge in v:
                w.add_edge(k, edge, [])

    print('saving mxtfile', mxtfile)
    with open(mxtfile, 'w') as fp:
        for k, v in pathy.hashdict.items():
            fp.write("%d,%s\n" % (k, " ".join(map(str, v))))
        fp.close()

    if args.label:
        print('note: used/assigned %d labels total' % (len(set(all_labels)),))
        print('counts:', label_counts)


if __name__ == '__main__':
    main()
