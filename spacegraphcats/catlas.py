#!/usr/bin/python
import sys, os, argparse
from os import path
from collections import deque, defaultdict
from graph import Graph
from rdomset import rdomset, calc_domination_graph, calc_dominators
from minhash import MinHash
from parser import parse_minhash, Writer

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]

class CAtlasBuilder:
    def __init__(self, graph, vsizes, domination, minhashes):
        self.graph = graph
        self.domgraph = domination.domgraph
        self.domset = domination.domset
        self.minhashes = minhashes

        self.assignment = domination.assignment
        self.sizes = vsizes

        self.level_threshold = 10
        self.minhash_size = 1000

    def set_level_threshold(self, thres):
        self.level_threshold = thres
        return self

    def set_minhash_size(self, size):
        self.minhash_size = size
        return self        

    def _build_component_catlas(self, comp, vertices):
        comp = set(comp)
          # Build first level 
        curr_level = {}
        curr_domgraph = self.domgraph.subgraph(comp)

        # Collect hashes of the assigned vertices 
        leaf_hashes = defaultdict(lambda: MinHash(self.minhash_size))        
        for u in vertices:
            for v in self.assignment[u]:
                leaf_hashes[v].merge(self.minhashes[u], self.minhash_size)

        # Create level 0
        for v in comp:
            curr_level[v] = CAtlas.create_leaf(v, self.sizes[v], leaf_hashes[v])

        # Build remaining levels
        level = 1 # Just for debugging
        while len(curr_level) > self.level_threshold:
            domset, augg = rdomset(curr_domgraph,1)
            dominators = calc_dominators(augg, domset, 1)
            next_domgraph, domset, dominators, next_assignment = calc_domination_graph(curr_domgraph, augg, domset, dominators, 1)

            child_map = defaultdict(set)
            for u in curr_domgraph:
                for v in next_assignment[u]:
                    child_map[v].add(u)

            next_level = {}
            # Debugging stats
            maxdeg, degsum, mindeg = 0, 0, float('inf')
            for v in next_domgraph:
                children = [curr_level[u] for u in child_map[v]]
                deg = len(children)
                maxdeg, mindeg = max(maxdeg, deg), min(mindeg, deg)
                degsum += deg
                next_level[v] = CAtlas.create_node(v, children, self.minhash_size)

            curr_domgraph = next_domgraph
            curr_level = next_level
            level += 1

        # Add root
        root = CAtlas.create_node('root', curr_level.values(), self.minhash_size)
        return root


    def build(self):
        comp_atlases = []

        # Compute data structure to handle components
        comp_lookup = self.domgraph.component_index() # Maps vertices to component ids
        components = defaultdict(set)
        for v in self.domgraph:
            components[comp_lookup[v]].add(v)
        components = list(components.values())

        num_comps = len(components)

        # Collect vertices from g that belong to the respective components
        print("Mapping graph vertices to dominators")
        graph_comps = defaultdict(set)
        n = len(self.graph)
        for i, u in enumerate(self.graph):
            if i % 10 == 0:
                print("\rProcessing node {}/{}".format(i,n), end="", flush=True)
            dominator = next(iter(self.assignment[u]))
            graph_comps[comp_lookup[dominator]].add(u)

        print("\rProcessing node {}/{}".format(n,n), end="", flush=True) # For my OCD

        print("\nBuilding catlasses for connected components")
        for i, comp in enumerate(components):
            print("\rProcessing component {}/{}".format(i,num_comps), end="", flush=True)
            comp_id = comp_lookup[next(iter(comp))] 
            comp_atlases.append(self._build_component_catlas(comp, graph_comps[comp_id]))

        if len(comp_atlases) == 1:
            return comp_atlases[0]

        curr_level = comp_atlases
        while len(curr_level) > self.level_threshold:
            next_level = []
            for block in chunks(curr_level, self.level_threshold):
                next_level.append(CAtlas.create_virtual_node(block, self.minhash_size))
            curr_level = next_level

        if len(curr_level) > 1:
            root = CAtlas.create_virtual_node(curr_level, self.minhash_size)
        else:
            root = curr_level[0]

        return root
      

class CAtlas:
    def __init__(self, id, level, size, children, minhash):
        self.id = id
        self.children = children
        self.size = size
        self.minhash = minhash
        self.level = level

    @staticmethod
    def create_leaf(id,  size, minhash):
        return CAtlas(id, 0, size, [], minhash)

    @staticmethod
    def create_node(id, children, hash_size):
        assert len(children) > 0
        size = 0
        level = next(iter(children)).level
        minhash = MinHash(hash_size)
        for c in children:
            assert c.level == level
            size += c.size
            minhash.merge(c.minhash,hash_size)
        return CAtlas(id, level+1, size, children, minhash)

    @staticmethod
    def create_virtual_node(children, hash_size):
        assert len(children) > 0
        size = 0
        maxlevel = max(children, key=lambda c: c.level).level
        minhash = MinHash(hash_size)
        for c in children:
            assert c.level <= maxlevel
            size += c.size
            minhash.merge(c.minhash,hash_size)

        return CAtlas('virtual', maxlevel+1, size, children, minhash)


    def score(self, Q, querysize):
        minsize = 1.0*min(len(Q),len(self.minhash))
        querysize = 1.0*len(Q)
        score = self.minhash.intersect(Q) / querysize
        return score

    def query(self, Q, querysize, threshold):
        res = set()
        score = self.score(Q,querysize)
        if score < threshold:
            print(score)
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

    def dfs(self): 
        yield self
        for c in self.children:
            for v in c.dfs():
                yield v

    def bfs(self): 
        curr_level = set([self])
        while len(curr_level) > 0:
            yield curr_level
            next_level = set()              
            for v in curr_level:
                next_level |= set(v.children)
            curr_level = next_level

    def write(self, projectpath, projectname, radius):
        """
            Write catlas to file
        """
        # Create unique ids for DAG-nodes
        idmap = {}
        id = 0
        for level in self.bfs():
            for v in level:
                idmap[v] = id
                id += 1

        # Store DAG in .gxt file
        with open(path.join(projectpath, "{}.catlas.{}.gxt".format(projectname, radius)), 'w') as f:
            writer = Writer(f, ['vertex', 'level'], [])
            for level in self.bfs():
                for v in level:
                    writer.add_vertex(idmap[v], v.size, [v.id, v.level])

            for level in self.bfs():
                for v in level:
                    for c in v.children:
                        writer.add_edge(idmap[v], idmap[c])
            writer.done()

        # Store hashes separately in .mxt
        with open(path.join(projectpath, "{}.catlas.{}.mxt".format(projectname, radius)), 'w') as f:
            for level in self.bfs():
                for v in level:
                    f.write(str(idmap[v])+',')
                    f.write(' '.join(map(str,v.minhash)))
                    f.write('\n')


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

