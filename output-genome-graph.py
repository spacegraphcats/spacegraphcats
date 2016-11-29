#! /usr/bin/env python3
import screed
import argparse
from sourmash_lib import MinHash
import os, os.path
from collections import OrderedDict, defaultdict
from spacegraphcats import graph_parser


KSIZE=31
CHUNKSIZE=int(1e5)

# minhash settings
MH_SIZE_DIVISOR=50
MH_MIN_SIZE=5

class Pathfinder(object):
    "Track segment IDs, adjacency lists, and MinHashes"
    def __init__(self, ksize, mxtfile, node_offset=0):
        self.ksize = ksize

        self.node_counter = 1 + node_offset
        self.nodes = {}                      # node IDs (int) to size
        self.nodes_to_kmers = {}             # node IDs (int) to kmers
        self.kmers_to_nodes = {}             # kmers to node IDs
        self.adjacencies = defaultdict(set)  # node to node
        self.labels = defaultdict(set)       # nodes to set of labels
        self.mxtfp = open(mxtfile, 'wt')

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
        self.nodes[node_id] = size

        kmer = min(visited)               # identify linear nodes by min(hash)
        self.nodes_to_kmers[node_id] = kmer
        self.kmers_to_nodes[kmer] = node_id

        return node_id

    def add_adjacency(self, node_id, adj):
        "Add an edge between two nodes to the cDBG."
        node_id, adj = min(node_id, adj), max(node_id, adj)
        
        x = self.adjacencies[node_id]
        x.add(adj)

    def add_label(self, kmer, label):
        x = self.labels[kmer]
        x.add(label)

    def add_minhash(self, path_id, mh):
        # save minhash info to disk
        mins = " ".join(map(str, mh.get_mins()))
        self.mxtfp.write('{0},{1}\n'.format(path_id, mins))


def chunks(seq, chunksize):
    for start in range(0, len(seq) - 1, chunksize):
        yield start, start+chunksize, seq[start:start+chunksize]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('seqfiles', nargs='+')
    p.add_argument('-o', '--output', default=None)
    args = p.parse_args()

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

    pathy = Pathfinder(KSIZE, mxtfile, 0)

    for seqfile in args.seqfiles:
        for record in screed.open(seqfile):
            for start, stop, chunk in chunks(record.sequence, 10000):
                node_id = pathy.new_linear_node([0], CHUNKSIZE)
                if start > 0:
                    pathy.add_adjacency(node_id - 1, node_id)

                mh = MinHash(CHUNKSIZE // MH_SIZE_DIVISOR, MH_MIN_SIZE)
                mh.add_sequence(chunk)
                pathy.add_minhash(node_id, mh)

    with open(gxtfile, 'w') as fp:
        w = graph_parser.Writer(fp, ['labels'], [])

        for k, v in pathy.nodes.items():
            w.add_vertex(k, v)

        for k, v in pathy.adjacencies.items():
            for edge in v:
                w.add_edge(k, edge, [])


if __name__ == '__main__':
    main()
