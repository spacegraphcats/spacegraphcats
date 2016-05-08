#! /usr/bin/env python
from __future__ import print_function
import khmer
import screed
import argparse
from collections import OrderedDict

DEFAULT_KSIZE=31
NODEGRAPH_SIZE=8e8 # small, big is 2e8

MH_SIZE_DIVISOR=50
MH_MIN_SIZE=5
MH_PRIME=9999999967

def kmers(sequence, K):
    for i in range(len(sequence) - K + 1):
        yield sequence[i:i+K]


def neighbors(kmer, cg):
    assert len(kmer) == cg.ksize()

    prefix = kmer[1:]
    postfix = kmer[:-1]

    for ch in 'ACGT':
        kmer = prefix + ch
        if cg.get(kmer):
            yield kmer

        kmer = ch + postfix
        if cg.get(kmer):
            yield kmer


class Pathfinder(object):
    def __init__(self, cg, bf):
        self.cg = cg
        self.bf = bf

        self.segment_counter = 1
        self.segments = {}
        self.segments_r = {}
        self.segment_start = {}
        self.adjacencies = {}
        self.adj_tips = {}
        self.hashdict = OrderedDict()

    def new_segment(self, kmer):
        if kmer in self.segments_r:
            return self.segments_r[kmer]

        this_id = self.segment_counter
        self.segment_counter += 1

        self.segments[this_id] = self.cg.ksize()
        self.segments_r[kmer] = this_id

        return this_id

    def new_linear_segment(self, size):
        this_id = self.segment_counter
        self.segment_counter += 1
        self.segments[this_id] = size
        return this_id

    def add_adjacency(self, node_id, adj):
        if adj < node_id:
            return
        
        x = self.adjacencies.get(node_id, set())
        x.add(adj)
        self.adjacencies[node_id] = x

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seqfiles', nargs='+')
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('-k', '--ksize', default=DEFAULT_KSIZE, type=int)
    parser.add_argument('-x', '--tablesize', default=NODEGRAPH_SIZE,
                            type=float)
    parser.add_argument('--force', action='store_true')
    parser.add_argument('--gml', action='store_true')
    args = parser.parse_args()

    assert args.ksize % 2, "ksize must be odd"
    assert args.output, "you probably want an output file"

    print('building graphs and loading files')
    cg = khmer.Countgraph(args.ksize, args.tablesize, 2)
    bf = khmer.Nodegraph(args.ksize, args.tablesize, 2)
    bf2 = khmer.Nodegraph(args.ksize, args.tablesize, 2)
    n = 0
    for seqfile in args.seqfiles:
        for record in screed.open(seqfile):
            n += 1
            if n % 10000 == 0:
                print('...', seqfile, n)
            cg.consume(record.sequence)

    fp_rate = khmer.calc_expected_collisions(cg, args.force, max_false_pos=.05)
    pathy = Pathfinder(cg, bf)

    print('finding high degree nodes')
    n = 0
    for seqfile in args.seqfiles:
        for record in screed.open(seqfile):
            n += 1
            if n % 10000 == 0:
                print('...2', seqfile, n)
            # walk across sequences, find all high degree nodes,
            # name them and cherish them.
            if min(bf2.get_kmer_counts(record.sequence)) == 0:
                bf2.consume(record.sequence)
                cg.find_high_degree_nodes(record.sequence)

    del bf2

    degree_nodes = cg.get_high_degree_nodes()
    for node in degree_nodes:
        pathy.new_segment(node)

    print('traversing linear segments from', len(degree_nodes), 'nodes')
    for n, k in enumerate(degree_nodes):
        if n % 10000 == 0:
            print('...', n, 'of', len(degree_nodes))
        k_id = pathy.segments_r[k]
        
        k_str = khmer.reverse_hash(k, cg.ksize())
        hashval = khmer._minhash.hash_murmur32(k_str) % MH_PRIME;
        pathy.hashdict[k_id] = [hashval]

        nbh = cg.neighbors(k)
        for nk in nbh:
            if cg.is_high_degree_node(nk):
                nk_id = pathy.segments_r[nk]
                pathy.add_adjacency(k_id, nk_id)
            else:
                size, conns, visited = cg.traverse(nk, bf)
                if size:
                    path_id = pathy.new_linear_segment(size)

                    v = [ khmer.reverse_hash(i, cg.ksize()) for i in visited ]
                    mh_size = max(len(visited) // MH_SIZE_DIVISOR, MH_MIN_SIZE)
                    mh = khmer.MinHash(mh_size, cg.ksize())
                    for kmer in v:
                        mh.add_sequence(kmer)
                    pathy.hashdict[path_id] = mh.get_mins()

                    for conn in conns:
                        conn_id = pathy.segments_r.get(conn)
                        pathy.add_adjacency(path_id, conn_id)
                        pathy.add_adjacency(conn_id, path_id)

    print(len(pathy.segments), 'segments, containing',
              sum(pathy.segments.values()), 'nodes')
    #for k in pathy.segments:
    #    print(k, pathy.segments[k], pathy.adjacencies.get(k, []))

    # save to GML
    if args.output:
        import parser

        print('saving to', args.output)
        fp = open(args.output, 'w')
        if args.gml:
            w = parser.GmlWriter(fp, [], [])
        else:
            w = parser.Writer(fp, [], [])
        
        for k, v in pathy.segments.items():
            w.add_vertex(k, v, [])

        for k, v in pathy.adjacencies.items():
            for edge in v:
                w.add_edge(k, edge, [])

        if not args.gml:
            fp = open(args.output + '.mxt', 'w')
            for k, v in pathy.hashdict.items():
                fp.write("%d,%s\n" % (k, " ".join(map(str, v))))
            fp.close()


if __name__ == '__main__':
    main()
