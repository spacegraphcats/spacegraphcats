#!/usr/bin/python
import sys
import os
import argparse
import networkx as nx
from collections import deque, defaultdict
from graph import Graph
from rdomset import rdomset, calc_domset_graph, calc_dominators
from minhash import MinHash
from parser import parse_minhash

class CAtlasBuilder:
    def __init__(self, G, hashes, hash_size, bottom_dist, upper_dist, top_level_size):
        self.levels = deque()
        self.bottom_dist = bottom_dist
        self.upper_dist = upper_dist
        self.top_level_size = top_level_size
        self.hashes = hashes
        self.hash_size = hash_size
        self.compute_levels(G)

    def compute_levels(self, G):
        """
            Build up the levels of the catlas from the bottom up
        """
        # get bottom level
        prev = len(G)
        H = self.create_new_level(G,bottom=True)
        # keep going until the top level is small
        while len(H) > self.top_level_size and len(H) < prev:
            prev = len(H)
            H = self.create_new_level(H)
        # set a root to point to the whole top level
        U = self.top_level_vertices()
        domset = {U[0]: U}
        if len(H) >= prev:
            self.levels[0] = domset
        else:
            self.levels.appendleft(domset)
        self.root = U[0]

    def create_new_level(self, G, bottom=False):
        """
            Create a new level and return the "minor" made from the
            intersection of dominated sets.
        """
        # use a different distance d on the bottom set
        if bottom:
            d = self.bottom_dist
        else:
            d = self.upper_dist
        # find the dominating set
        domset, G_aug = rdomset(G,d)
        dominators = calc_dominators(G_aug, domset, d)
        level = defaultdict(list)
        # make lookup of dominated by dominators
        for v, dom_by_radius in dominators.iteritems():
            for r, S in dom_by_radius.iteritems():
                for u in S:
                    level[u].append(v)
        # add structure
        self.levels.appendleft(level)
        return calc_domset_graph(domset, dominators, d)

    def traverse(self,cutoff=None):
        """
            Traverse the catlas up until cutoff level
        """
        # cutoff at the bottom if not otherwise specified
        if cutoff is None:
            cutoff = self.depth()-1
        i = 0
        # continue down levels until we reach the cutoff
        while (i <= cutoff):
            print "Level", i, len(self.levels[i]), max(len(S) for _,S in self.levels[i].iteritems())
            i += 1

    def top_level_vertices(self):
        """
            Return the vertices at the top level
        """
        return self.levels[0].keys()

    def depth(self):
        """
            return the depth of the catlas
        """
        return len(self.levels)

    def build_minhash(self,seq,base=False):
        """
            Build a minhash from a sequence of minhashes
            in seq
        """
        h = None
        if base:
            h = MinHash(self.hash_size)
            for v in seq:
                h.add(self.hashes[v])
        else:
            h = reduce( lambda h1,h2: h1.merge(h2,self.hash_size), [c.minhash for c in seq] )
            assert len(h) <= self.hash_size
        return h

    def build_catlas(self):
        return self.recursive_build(self.root,0)

    def recursive_build(self, id, level):
        """
            Build the CAtlas object recursively
        """
        # base case of recursion is first domset
        if level == self.depth()-1:
            children = []
            minhash = self.build_minhash(self.levels[level][id], base=True)
            size = len(self.levels[level][id])
        else:
            # recursively build children
            children = [self.recursive_build(i, level+1) for i in self.levels[level][id]]
            minhash = self.build_minhash(children)
            size = sum(c.size for c in children)
        return CAtlas(id, level, size, children, minhash)

    @classmethod
    def from_file(cls, input):
        """
            build an catlas from a file dump
        """
        #TODO:  make sure this is a valid catlas tree
        raw_data = {}
        # extract all the data on the first pass
        for line in input:
            line = line.strip()
            if len(line) == 0:
                continue
            id, level, size_str, children_str, mh_str  = line.split(';')
            raw_data[(int(level),int(id))] = (size_str, children_str, mh_str)

        # find root
        _, root = min(raw_data.keys())
        # build recursively
        return cls.recursive_build_from_file(raw_data, root, 0, 0)

    @classmethod
    def recursive_build_from_file(cls, raw_data, id, level, hash_size):
        """
            build catlases recursively from the file information
        """
        size_str, children_str, mh_str = raw_data[(level, id)]
        size = int(size_str)
        children = [] if len(children_str) == 0 else map(int, children_str.split(','))

        children = [cls.recursive_build_from_file(raw_data, c, level+1, hash_size) \
                    for c in children]
        mh_list = [int(h) for h in mh_str.split(',')]

        minhash = MinHash.from_list(mh_list)

        return CAtlas(id, level, size, children, minhash)


class CAtlas:
    def __init__(self, id, level, size, children, minhash):
        self.id = id
        self.children = children
        self.size = size
        self.minhash = minhash
        self.level = level

    def score(self, Q, querysize):
        minsize = 1.0*min(len(Q),len(self.minhash))
        querysize = 1.0*len(Q)
        score = len(self.minhash.intersect(Q)) / querysize
        return score

    def query(self, Q, querysize, threshold):
        res = set()
        score = self.score(Q,querysize)
        if score < threshold:
            print score
            return res

        if len(self.children) == 0:
            res.add(self.id)
            return res

        for c in self.children:
            res |= c.query(Q, querysize, threshold)

        return res

    def leaves(self):
        if len(self.children) == 0:
            return set([self.id])
        res = set()
        for c in self.children:
            res |= c.leaves()
        return res

    def leaves2(self):
        if len(self.children) == 0:
            return set([self])
        res = set()
        for c in self.children:
            res |= c.leaves2()
        return res


    def to_file(self, stream):
        """
            write properties to filestream stream
        """
        children_id_str = [str(c.id) for c in self.children]
        children_str = ','.join(children_id_str)
        stream.write("{};{};{};{};{}\n".format(self.id, self.level, self.size, children_str, self.minhash))
        for c in self.children:
            c.to_file(stream)

def read_minhashes(inputfile):
    res = {}

    def add_hash(v, hashes):
        res[int(v)] = MinHash.from_list(hashes)

    parse_minhash(inputfile, add_hash)
    return res

def main(inputgraph, inputhash, output, hash_size, d_bottom, d_upper, top_level_size, test):
    G,node_attrs,edge_attrs = Graph.from_gxt(inputgraph)
    hashes = read_minhashes(inputhash)

    for v in G:
        assert v in hashes, "{} has no minhash".format(v)

    G.remove_loops()
    AB = CAtlasBuilder(G, hashes, hash_size, d_bottom, d_upper, top_level_size)
    print("Input graph has size", len(G))
    catlas = AB.build_catlas()
    catlas.to_file(output)

    print("CAtlas has {} leaves".format(len(catlas.leaves())))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputgraph', help='gxt input file', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('inputhash', help='mxt input file', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('output', help='catlas output file', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument("hashsize", help="number of hashes to keep", type=int)
    parser.add_argument("--bottom", help="bottom level domset distance", type=int, default=10)
    parser.add_argument("--inner", help="upper levels domset distance", type=int, default=1)
    parser.add_argument("--t", help="top level size threshold", type=int, default=5)
    parser.add_argument('--test', action='store_true')
    args = parser.parse_args()

    main(args.inputgraph, args.inputhash, args.output, args.hashsize, args.bottom, args.inner, args.t, args.test)

