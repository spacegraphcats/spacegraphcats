from __future__ import print_function

import sys
import khmer
from sourmash_lib import MinHash, MAX_HASH
import screed
from collections import OrderedDict, defaultdict
import os, os.path
from spacegraphcats import graph_parser


# minhash settings
SCALED=100.0

class Pathfinder(object):
    "Track segment IDs, adjacency lists, and MinHashes"
    def __init__(self, ksize, gxtfile, mxtfile):
        self.ksize = ksize

        self.node_counter = 0
        self.nodes = {}                      # node IDs (int) to size
        self.nodes_to_kmers = {}             # node IDs (int) to kmers
        self.kmers_to_nodes = {}             # kmers to node IDs
        self.adjacencies = defaultdict(set)  # node to node
        self.labels = defaultdict(set)       # nodes to set of labels
        self.mxtfp = open(mxtfile, 'wt')
        self.adjfp = open(gxtfile + '.adj', 'wt')

    def new_hdn(self, kmer):
        "Add a new high-degree node to the cDBG."
        if kmer in self.kmers_to_nodes:
            return self.kmers_to_nodes[kmer]

        this_id = self.node_counter
        self.node_counter += 1

        self.nodes[this_id] = 1
        self.nodes_to_kmers[this_id] = kmer
        self.kmers_to_nodes[kmer] = this_id

        return this_id

    def new_linear_node(self, visited, size):
        "Add a new linear path to the cDBG."
        node_id = self.node_counter
        self.node_counter += 1

        return node_id

    def add_adjacency(self, node_id, adj):
        "Add an edge between two nodes to the cDBG."
        node_id, adj = min(node_id, adj), max(node_id, adj)

        self.adjfp.write('{},{}\n'.format(node_id, adj))
        
    def add_label(self, kmer, label):
        x = self.labels[kmer]
        x.add(label)

    def add_minhash(self, path_id, mh):
        # save minhash info to disk
        mins = " ".join(map(str, mh.get_mins()))
        self.mxtfp.write('{0} {1}\n'.format(path_id, mins))


def traverse_and_mark_linear_paths(graph, nk, stop_bf, pathy, degree_nodes):
    size, adj_kmers, visited = graph.traverse_linear_path(nk, degree_nodes,
                                                          stop_bf)
    if not size:                          # 0 length paths
        return

    # get an ID for the new path
    path_id = pathy.new_linear_node(visited, size)

    # add all adjacencies
    for kmer in adj_kmers:
        adj_node_id = pathy.kmers_to_nodes[kmer]
        pathy.add_adjacency(path_id, adj_node_id)

    # next, calculate minhash from visited k-mers
    v = [ khmer.reverse_hash(i, graph.ksize()) for i in visited ]

    mh = MinHash(n=0, ksize=graph.ksize(), max_hash=int(MAX_HASH / SCALED))
    for kmer in v:
        mh.add_sequence(kmer)

    print('XXX', len(v), len(mh.get_mins()))

    pathy.add_minhash(path_id, mh)


def run(args):

    # @CTB this is kind of a hack - nothing tricky going on, just want to
    # specify memory on the command line rather than graph size...
    graph_tablesize = int(args.memory * 8.0 / 4.0)

    assert args.ksize % 2, "ksize must be odd"

    if args.label:
        label_list = []

    output_dir = args.output
    if not output_dir:
        if len(args.seqfiles) > 1:
            print('** please specify an output directory with -o',
                  file=sys.stderr)
            sys.exit(-1)

        output_dir = os.path.basename(args.seqfiles[0])
        if output_dir.endswith('.fa'):
            output_dir = output_dir[:-3]
        elif output_dir.endswith('.fa.gz'):
            output_dir = output_dir[:-6]

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
    if args.loadgraph:
        print('loading nodegraph from:', args.loadgraph)
        graph = khmer.load_nodegraph(args.loadgraph)
        print('creating accompanying stopgraph')
        ksize = graph.ksize()
        hashsizes = graph.hashsizes()
        stop_bf = khmer._Nodegraph(ksize, hashsizes)
    else:
        print('building graphs and loading files')

        # Create graph and a stop bloom filter - one for loading, one for
        # traversing. Create them all here so that we can error out quickly
        # if memory is a problem.

        graph = khmer.Nodegraph(args.ksize, graph_tablesize, 2)
        stop_bf = khmer.Nodegraph(args.ksize, graph_tablesize, 2)
        n = 0

        # load in all of the input sequences, one file at a time.
        for seqfile in args.seqfiles:
            for record in screed.open(seqfile):
                if len(record.sequence) < graph.ksize(): continue
                n += 1
                if n % 10000 == 0:
                    print('...', seqfile, n)
                graph.consume(record.sequence)

        # complain if too small set of graphs was used.
        fp_rate = khmer.calc_expected_collisions(graph,
                                                 args.force, max_false_pos=.05)

    ksize = graph.ksize()

    # initialize the object that will track information for us.
    pathy = Pathfinder(ksize, gxtfile, mxtfile)

    print('finding high degree nodes')
    if args.label:
        print('(and labeling them, per request)')
    degree_nodes = khmer.HashSet(ksize)
    n = 0
    for seqfile in args.seqfiles:
        for record in screed.open(seqfile):
            if len(record.sequence) < ksize: continue
            n += 1
            if n % 10000 == 0:
                print('...2', seqfile, n)
            # walk across sequences, find all high degree nodes,
            # name them and cherish them.
            these_hdn = graph.find_high_degree_nodes(record.sequence)
            degree_nodes += these_hdn
            if args.label:
                label_list.append(record.name)
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
        k_str = khmer.reverse_hash(k, ksize)
        mh = MinHash(n=0, ksize=ksize, max_hash=round(MAX_HASH / SCALED))
        mh.add_sequence(k_str)
        pathy.add_minhash(k_id, mh)

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

    del graph
    del stop_bf

    # save to GXT/MXT.
    print('saving gxtfile', gxtfile)

    all_labels = set()
    label_counts = {}

    pathy.adjfp.close()
    adj_filename = open(gxtfile + '.adj', 'rt')

    # this uniqifies the edges...
    for line in adj_filename:
        a, b = line.split(',')
        a = int(a)
        b = int(b)
        pathy.adjacencies[a].add(b)

    try:
        os.unlink(gxtfile + '.adj')
    except:
        print('cannot remove', gxtfile + '.adj')

    # ...and now print them out.
    edges = []
    for k, v in pathy.adjacencies.items():
        for dest in v:
            # don't add loops
            if (k != dest):
                edges.append((k, dest))

    graph_parser.write(open(gxtfile, 'wt'), pathy.node_counter, edges)

    if args.label:
        print('note: used/assigned %d labels total' % (len(set(all_labels)),))
        print('counts:', label_counts)

        assert label_list
        print('dumping label list now.')
        label_file = os.path.basename(output_dir) + '.labels.txt'
        label_file = os.path.join(output_dir, label_file)

        with open(label_file, "wt") as fp:
            for n, label in enumerate(label_list):
                fp.write("{} {}\n".format(n + 0, label))
