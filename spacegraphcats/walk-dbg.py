#! /usr/bin/env python3

import khmer
import screed
import argparse
from collections import OrderedDict

# graph settings
DEFAULT_KSIZE=31
NODEGRAPH_SIZE=8e8 # small, big is 2e8

# minhash settings
MH_SIZE_DIVISOR=50
MH_MIN_SIZE=5
MH_PRIME=9999999967

class Pathfinder(object):
    "Track segment IDs, adjacency lists, and MinHashes"
    def __init__(self, ksize):
        self.ksize = ksize

        self.segment_counter = 1
        self.segments = {}                # segment IDs (int) to size
        self.segments_f = {}              # segment IDs (int) to kmers
        self.segments_r = {}              # kmers to segment IDs
        self.adjacencies = {}
        self.hashdict = OrderedDict()
        self.labels = {}

    def new_segment(self, kmer):
        if kmer in self.segments_r:
            return self.segments_r[kmer]

        this_id = self.segment_counter
        self.segment_counter += 1

        self.segments[this_id] = self.ksize
        self.segments_f[this_id] = kmer
        self.segments_r[kmer] = this_id

        return this_id

    def new_linear_segment(self, size):
        this_id = self.segment_counter
        self.segment_counter += 1
        self.segments[this_id] = size
        return this_id

    def add_adjacency(self, node_id, adj):
        node_id, adj = min(node_id, adj), max(node_id, adj)
        
        x = self.adjacencies.get(node_id, set())
        x.add(adj)
        self.adjacencies[node_id] = x

    def add_label(self, node_id, label):
        x = self.labels.get(node_id, [])
        x.append(label)
        self.labels[node_id] = x


def traverse_and_mark_linear_paths(graph, nk, stop_bf, pathy, degree_nodes):
    size, conns, visited = graph.traverse_linear_path(nk, degree_nodes,
                                                      stop_bf)
    if not size:
        return

    # give it a segment ID
    path_id = pathy.new_linear_segment(size)

    # for all adjacencies, add.
    for conn in conns:
        conn_id = pathy.segments_r.get(conn)
        pathy.add_adjacency(path_id, conn_id)

    # next, calculate minhash from visited k-mers
    v = [ khmer.reverse_hash(i, graph.ksize()) for i in visited ]
    mh_size = max(len(visited) // MH_SIZE_DIVISOR, MH_MIN_SIZE)

    mh = khmer.MinHash(mh_size, graph.ksize())
    for kmer in v:
        mh.add_sequence(kmer)

    # save minhash info
    pathy.hashdict[path_id] = mh.get_mins()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqfiles', nargs='+')
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('-k', '--ksize', default=DEFAULT_KSIZE, type=int)
    parser.add_argument('-x', '--tablesize', default=NODEGRAPH_SIZE,
                            type=float)
    parser.add_argument('--force', action='store_true')
    parser.add_argument('--gml', action='store_true')
    parser.add_argument('--label', action='store_true')
    args = parser.parse_args()

    assert args.ksize % 2, "ksize must be odd"
    assert args.output, "you probably want an output file"

    print('building graphs and loading files')

    # Create graph, and two stop bloom filters - one for loading, one for
    # traversing. Create them all here so that we can error out quickly
    # if memory is a problem.

    graph = khmer.Nodegraph(args.ksize, args.tablesize, 2)
    stop_bf = khmer.Nodegraph(args.ksize, args.tablesize, 2)
    stop_bf2 = khmer.Nodegraph(args.ksize, args.tablesize, 2)
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
    pathy = Pathfinder(args.ksize)

    print('finding high degree nodes')
    degree_nodes = khmer.HashSet(args.ksize)
    n = 0
    for seqfile in args.seqfiles:
        for record in screed.open(seqfile):
            if len(record.sequence) < args.ksize: continue
            n += 1
            if n % 10000 == 0:
                print('...2', seqfile, n)
            # walk across sequences, find all high degree nodes,
            # name them and cherish them. Don't do this on identical sequences.
            if min(stop_bf2.get_kmer_counts(record.sequence)) == 0:
                stop_bf2.consume(record.sequence)
                degree_nodes += graph.find_high_degree_nodes(record.sequence)
                if args.label:
                    for node_id in degree_nodes:
                        pathy.add_label(node_id, n)
    del stop_bf2

    # get all of the degree > 2 nodes and give them IDs.
    for node in degree_nodes:
        pathy.new_segment(node)

    print('traversing linear segments from', len(degree_nodes), 'nodes')

    # now traverse from each high degree nodes into all neighboring nodes,
    # seeking adjacencies.  if neighbor is high degree node, add it to
    # adjacencies; if neighbor is not, then traverse the linear path.  also
    # track minhashes while we're at it.
    for n, k in enumerate(degree_nodes):
        if n % 10000 == 0:
            print('...', n, 'of', len(degree_nodes))

        # retrieve the segment ID of the primary node.
        k_id = pathy.segments_r[k]

        # add its hash value.
        k_str = khmer.reverse_hash(k, graph.ksize())
        hashval = khmer._minhash.hash_murmur32(k_str) % MH_PRIME;
        pathy.hashdict[k_id] = [hashval]

        # find all the neighbors of this high-degree node.
        nbh = graph.neighbors(k)
        for nk in nbh:
            # neighbor is high degree? fine, mark its adjacencies.
            if nk in degree_nodes:
                nk_id = pathy.segments_r[nk]
                pathy.add_adjacency(k_id, nk_id)
            else:
                # linear! walk it.
                traverse_and_mark_linear_paths(graph, nk, stop_bf, pathy,
                                               degree_nodes)

    print(len(pathy.segments), 'segments, containing',
              sum(pathy.segments.values()), 'nodes')

    # save to GXT or GML.
    if args.output:
        import parser

        print('saving graph to', args.output)
        fp = open(args.output, 'w')
        if args.gml:
            w = parser.GmlWriter(fp, [], [])
        else:
            w = parser.Writer(fp, ['labels'], [])
        
        for k, v in pathy.segments.items():
            kmer = pathy.segments_f.get(k)
            l = ""
            if kmer:
                l = pathy.labels.get(kmer, "")
                if l:
                    l = " ".join(map(str, l))
            w.add_vertex(k, v, [l])

        for k, v in pathy.adjacencies.items():
            for edge in v:
                w.add_edge(k, edge, [])

        if not args.gml:
            parts = args.output.split('.')
            if parts[-1] == 'gxt':
                parts[-1] = 'mxt'
            else:
                parts.append('mxt')
            fp = open('.'.join(parts), 'w')
            print('saving minhashes to', ".".join(parts))
            for k, v in pathy.hashdict.items():
                fp.write("%d,%s\n" % (k, " ".join(map(str, v))))
            fp.close()


if __name__ == '__main__':
    main()
