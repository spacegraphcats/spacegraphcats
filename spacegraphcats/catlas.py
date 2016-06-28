#! /usr/bin/env python3
from __future__ import print_function

import sys, os, argparse
from os import path
from collections import deque, defaultdict

from spacegraphcats.graph import Graph, VertexDict
from spacegraphcats.rdomset import rdomset, calc_domination_graph, calc_dominators
from spacegraphcats.minhash import MinHash
from spacegraphcats.graph_parser import parse_minhash, Writer

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
                print("\rProcessing node {}/{}".format(i,n), end="")
                sys.stdout.flush()
            dominator = next(iter(self.assignment[u]))
            graph_comps[comp_lookup[dominator]].add(u)

        print("\rProcessing node {}/{}".format(n,n), end="") # For my OCD
        sys.stdout.flush()

        print("\nBuilding catlasses for connected components")
        for i, comp in enumerate(components):
            print("\rProcessing component {}/{}".format(i,num_comps), end="")
            sys.stdout.flush()
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
    class Scoring:
        """ Given a search frontier, a minhash-query and a threshold,
            return a score (higher == better) of how well the frontier
            represents the query.
        """ 
        @staticmethod
        def jaccard(frontier, mhquery, threshold):
            covered = set()
            for node in frontier:
                covered |= set(node.minhash)
            
            inter = len(set(mhquery) & covered)
            union = len(set(mhquery) | covered)

            return inter / union

        @staticmethod
        def jaccard_leaves(frontier, mhquery, threshold):
            leaves = set()
            for node in frontier:
                leaves |= node.leaves2()

            covered = set()
            for node in leaves:
                covered |= set(node.minhash)
            
            inter = len(set(mhquery) & covered)
            union = len(set(mhquery) | covered)

            return inter / union            

        @staticmethod
        def height(frontier, mhquery, threshold):
            max_level = max(frontier, key=lambda node: node.level).level
            return 1/(max_level+1)

        @staticmethod
        def avg_height(frontier, mhquery, threshold):
            height_sum = sum([node.level for node in frontier])
            avg_height = height_sum / len(frontier)
            return 1/(avg_height+1)

        @staticmethod
        def height_jaccard(frontier, mhquery, threshold):
            return CAtlas.Scoring.jaccard(frontier, mhquery, threshold) \
                  * CAtlas.Scoring.avg_height(frontier, mhquery, threshold)

    class Selection:
        """ Given a search frontier, a minhash-query and a threshold,
            return the 'worst' node of the frontier, that is, the node
            that will most likely benefit from refinement.
        """
        @staticmethod
        def smallest_intersection(frontier, mhquery, threshold):
            # candidates = list(filter(lambda node: node.level != 0, frontier))
            # if len(candidates) != 0:
            #     return min(candidates, key=lambda node: mhquery.intersect(node.minhash))
            # else:
            return min(frontier, key=lambda node: mhquery.intersect(node.minhash))

        @staticmethod
        def highest_smallest_intersection(frontier, mhquery, threshold):
            max_level = max(frontier, key=lambda node: node.level).level
            candidates = list(filter(lambda node: node.level == max_level, frontier))
            return CAtlas.Selection.smallest_intersection(candidates, mhquery, threshold)

    class Refinement:
        """ Given a 'bad node' and a minhash-query and a threshold,
            return a subset of that node's children that covers the query as well as possible.
        """            
        @staticmethod
        def greedy_coverage(bad_node, mhquery, threshold):
            uncovered = set()
            for c in bad_node.children:
                uncovered |= set(c.minhash)
            uncovered &= set(mhquery)

            candidates = set(bad_node.children)
            res = set()
            while len(uncovered) > 0:
                # Greedily pick child with largest intersection of uncovered hashes
                best = max(candidates, key=lambda node: len(uncovered & set(node.minhash)) )
                candidates.remove(best)
                res.add(best)
                uncovered -= set(best.minhash)
            return res

        @staticmethod
        def greedy_coverage_leaves(bad_node, mhquery, threshold):
            leaves = bad_node.leaves2()

            uncovered = set()
            for c in leaves:
                uncovered |= set(c.minhash)
            uncovered &= set(mhquery)

            candidates = set(bad_node.children)
            res = set()
            while len(uncovered) > 0 and len(candidates) > 0:
                # Greedily pick child with largest intersection of uncovered hashes
                best = max(candidates, key=lambda node: len(uncovered & set(node.minhash)) )
                candidates.remove(best)
                res.add(best)
                uncovered -= set(best.minhash)
            return res            
            

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

    def query(self, mhquery, threshold, scoring_strat, selection_strat, refinement_strat):
        frontier = set([self])

        while True:
            print("Current frontier", [(n.id, n.level) for n in frontier])            
            scoreBefore = scoring_strat(frontier, mhquery, threshold)
            print("Current score {:.5f}".format(scoreBefore))

            bad_node = selection_strat(frontier, mhquery, threshold)
            print("Refining node {} with {} children".format(bad_node.id, len(bad_node.children)))

            if bad_node.level == 0:
                frontier.remove(bad_node)
                refinement = set()
            else:
                frontier.remove(bad_node)
                refinement = set(refinement_strat(bad_node, mhquery, threshold)) 
                frontier |= refinement                
                print("Choose {} out of {} children for refinement".format(len(refinement), len(bad_node.children)))

            scoreAfter = scoring_strat(frontier, mhquery, threshold)            
            print("Score after refinement {:.5f}".format(scoreAfter))

            if scoreBefore > scoreAfter:
                # Improvement step failed: roll back and stop.
                frontier -= refinement
                frontier.add(bad_node)
                if bad_node.level == 0:
                    print("Stopped at bad leaf.")
                break

        res = set()
        for node in frontier:
            res |= node.leaves()

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

    def write(self, projectpath, projectname, radius, offset=0):
        """
            Write catlas to file
        """
        # Create unique ids for DAG-nodes
        idmap = {}
        id = offset
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

    @staticmethod 
    def read(catlas_gxt, catlas_mxt, radius):
        graph, node_attr, _ = Graph.from_gxt(open(catlas_gxt, 'r'))
        hashes = VertexDict.from_mxt(open(catlas_mxt, 'r'))
        # hashes = defaultdict(int)

        # Parse attributes and partition catlas nodes into layers
        levels = defaultdict(set)
        for v in graph:
            node_attr[v]['level'] = int(node_attr[v]['level'])
            node_attr[v]['size'] = int(node_attr[v]['size'])
            l = node_attr[v]['level']
            levels[l].add(v)

        # Build catlas
        cat_nodes = {}
        for v in levels[0]:
            size = node_attr[v]['size']
            id = int(node_attr[v]['vertex'])
            cat_nodes[v] = CAtlas(id, 0, size, [], hashes[v])

        max_level = max(levels.keys())
        assert(len(levels[max_level]) == 1)
        for l in range(1, max_level+1):
            for v in levels[l]:
                size = node_attr[v]['size']
                children = filter( lambda x: node_attr[x]['level'] < l, graph.neighbours(v))
                children = list(map( lambda x: cat_nodes[x], children))
                id = node_attr[v]['vertex']
                try:
                    id = int(id)
                except ValueError:
                    pass
                cat_nodes[v] = CAtlas(id, l, size, children, hashes[v])

        root_id = next(iter(levels[max_level]))
        return cat_nodes[root_id]

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
    parser.add_argument('project', help='project directory', type=str )
    args = parser.parse_args()

    print(args.project)

