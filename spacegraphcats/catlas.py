#! /usr/bin/env python3
from __future__ import print_function

import sys, os, argparse
from os import path
from collections import deque, defaultdict
from operator import itemgetter

from spacegraphcats.graph import Graph, VertexDict
from spacegraphcats.rdomset import rdomset, calc_domination_graph, calc_dominators
from spacegraphcats.graph_parser import parse_minhash, Writer
from sourmash_lib import MinHash

KSIZE=31

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
        if minhashes:
            self.minhash_size = 1000
        else:
            self.minhash_size = 0

    def set_level_threshold(self, thres):
        self.level_threshold = thres
        return self

    def set_minhash_size(self, size):
        self.minhash_size = size
        return self        

    def _build_component_catlas(self, comp, vertices, start_id):
        comp = set(comp)
          # Build first level 
        curr_level = {}
        curr_domgraph = self.domgraph.subgraph(comp)

        # Collect hashes of the assigned vertices 
        leaf_hashes = defaultdict(lambda: MinHash(self.minhash_size, KSIZE))
        if self.minhashes:
            for u in vertices:
                for v in self.assignment[u]:
                    leaf_hashes[v].merge(self.minhashes[u])

        # Create level 0
        id = start_id
        for v in comp:
            curr_level[v] = CAtlas.create_leaf(id, v, self.sizes[v], leaf_hashes.get(v))
            id += 1

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
                next_level[v] = CAtlas.create_node(id, v, children, self.minhash_size)
                id += 1

            curr_domgraph = next_domgraph
            curr_level = next_level
            level += 1

        # Add root
        root = CAtlas.create_node(id, 'root', curr_level.values(), self.minhash_size)
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
        id = 0
        for i, comp in enumerate(components):
            print("\rProcessing component {}/{}".format(i,num_comps), end="")
            sys.stdout.flush()
            comp_id = comp_lookup[next(iter(comp))] 
            comp_atlases.append(self._build_component_catlas(comp, graph_comps[comp_id], id))
            id = comp_atlases[-1].id + 1

        if len(comp_atlases) == 1:
            return comp_atlases[0]

        curr_level = comp_atlases
        while len(curr_level) > self.level_threshold:
            next_level = []
            for block in chunks(curr_level, self.level_threshold):
                next_level.append(CAtlas.create_virtual_node(id, block, self.minhash_size))
                id += 1
            curr_level = next_level

        if len(curr_level) > 1:
            root = CAtlas.create_virtual_node(id, curr_level, self.minhash_size)
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
            """
                Returns the jaccard similarity between the query and the frontier
            """
            # make a minhash for the frontier in its entirety
            frontier_mh = MinHash(len(mhquery), KSIZE)
            for node in frontier:
                frontier_mh.merge(node.minhash)
            # compute jaccard similarity
            jacc = mhquery.compare(frontier_mh)
            return jacc

        @staticmethod
        def weighted_jaccard(frontier, mhquery, threshold):
            # not sure weighted jaccard is possible using minhashes
            raise NotImplementedError
            """
            frontier_mh = MinHash(len(mhquery), KSIZE)
            fuzz = 0
            for node in frontier:
                inters = node.minhash.count_common(mhquery)
                inter_perc = inters / node.size
                fuzz += (1-inter_perc) * node.size
                frontier_mh.merge(node.minhash)
            
            inter = len(qset & covered) 
            union = len(qset | covered) + fuzz

            return inter / union    
            """

        @staticmethod
        def jaccard_leaves(frontier, mhquery, threshold):
            leaves = set()
            for node in frontier:
                leaves |= node.leaves2()

            # make a minhash for the leaves 
            leaves_mh = MinHash(len(mhquery), KSIZE)
            for node in leaves:
                leaves_mh.merge(node.minhash)
            # compute jaccard similarity
            jacc = mhquery.compare(leaves_mh)
            return jacc         

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
        def worst_node(frontier, mhquery, threshold):
            """
            return the smallest jaccard similarity with the query
            """
            def jacc(x):
                return mhquery.compare(x.minhash)
            return min(map(jacc, frontier))

        @staticmethod
        def height_jaccard(frontier, mhquery, threshold):
            return CAtlas.Scoring.jaccard(frontier, mhquery, threshold) \
                  * CAtlas.Scoring.avg_height(frontier, mhquery, threshold)

        @staticmethod
        def height_weighted_jaccard(frontier, mhquery, threshold):
            return CAtlas.Scoring.weighted_jaccard(frontier, mhquery, threshold) \
                  * CAtlas.Scoring.avg_height(frontier, mhquery, threshold)**0.5

    class Selection:
        """ Given a search frontier, a minhash-query and a threshold,
            return the 'worst' node of the frontier, that is, the node
            that will most likely benefit from refinement.
        """
        @staticmethod
        def smallest_intersection(frontier, mhquery, threshold):
            return min(frontier, key=lambda node: mhquery.compare(node.minhash))

        @staticmethod
        def largest_weighted_intersection(frontier, mhquery, threshold):
            return max(frontier, key=lambda node: mhquery.compare(node.minhash)*node.size)

        @staticmethod
        def largest_weighted_intersection_cached():
            weighted_intersection_cache = {}
            def lwi(frontier, mhquery, threshold):
                largest_node = None
                largest_score = -1
                for node in frontier:
                    # compute comparison if it's not in cache already
                    if node not in weighted_intersection_cache:
                        val = mhquery.compare(node.minhash)*node.size
                        weighted_intersection_cache[node] = val
                    # check if largest
                    score = weighted_intersection_cache[node]
                    if score > largest_score:
                        largest_score = score
                        largest_node = node
                assert largest_node is not None
                return largest_node
            return lwi

        @staticmethod
        def largest_intersection_height(frontier, mhquery, threshold):
            return max(frontier, key=lambda node: mhquery.compare(node.minhash) * node.size * (1+node.level) )
        
        @staticmethod
        def weighted_jaccard(frontier, mhquery, threshold):
            def weighted_jacc(node):
                return mhquery.compare(node.minhash) * node.size             
            return min(frontier, key=weighted_jacc)     

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
                uncovered |= set(c.minhash.get_mins())
            uncovered &= set(mhquery.get_mins())

            candidates = set(filter(lambda n: mhquery.compare(n.minhash) >= threshold 
                                           or n.minhash.compare(mhquery) >= threshold, bad_node.children))
            res = set()
            while len(uncovered) > 0 and len(candidates) > 0:
                # Greedily pick child with largest intersection of uncovered hashes
                best = max(candidates, key=lambda node: len(uncovered & set(node.minhash.get_mins())) )
                candidates.remove(best)
                res.add(best)
                uncovered -= set(best.minhash.get_mins())
            return res

        @staticmethod
        def greedy_coverage_leaves(bad_node, mhquery, threshold):
            leaves = bad_node.leaves2()

            uncovered = set()
            for c in leaves:
                uncovered |= set(c.minhash.get_mins())
            uncovered &= set(mhquery.get_mins())

            candidates = set(bad_node.children)
            res = set()
            while len(uncovered) > 0 and len(candidates) > 0:
                # Greedily pick child with largest intersection of uncovered hashes
                best = max(candidates, key=lambda node: len(uncovered & set(node.minhash.get_mins())) )
                candidates.remove(best)
                res.add(best)
                uncovered -= set(best.minhash.get_mins())
            return res            
            

    # CAtlas.__init__
    def __init__(self, id, vertex, level, size, children, minhash):
        self.id = id
        self.vertex = vertex
        self.children = children
        self.size = size
        self.minhash = minhash
        self.level = level

    @staticmethod
    def create_leaf(id, vertex, size, minhash):
        return CAtlas(id, vertex, 0, size, [], minhash)

    @staticmethod
    def create_node(id, vertex, children, hash_size):
        assert len(children) > 0
        size = 0
        level = next(iter(children)).level
        if hash_size:
            minhash = MinHash(hash_size, KSIZE)
            for c in children:
                assert c.level == level
                size += c.size
                minhash.merge(c.minhash)
        else:
            minhash = None
        return CAtlas(id, vertex, level+1, size, children, minhash)

    @staticmethod
    def create_virtual_node(id, children, hash_size):
        assert len(children) > 0
        size = 0
        maxlevel = max(children, key=lambda c: c.level).level
        if hash_size:
            minhash = MinHash(hash_size, KSIZE)
            for c in children:
                assert c.level <= maxlevel
                size += c.size
                minhash.merge(c.minhash)
        else:
            minhash = None

        return CAtlas(id, 'virtual', maxlevel+1, size, children, minhash)

    def query(self, mhquery, threshold, scoring_strat, selection_strat, refinement_strat):
        frontier = set([self])

        while True:
            scoreBefore = scoring_strat(frontier, mhquery, threshold)
            bad_node = selection_strat(frontier, mhquery, threshold)

            if bad_node.level == 0:
                frontier.remove(bad_node)
                refinement = set()
            else:
                frontier.remove(bad_node)
                refinement = set(refinement_strat(bad_node, mhquery, threshold)) 
                frontier |= refinement                

            scoreAfter = scoring_strat(frontier, mhquery, threshold)            

            if scoreBefore > scoreAfter:
                frontier -= refinement
                frontier.add(bad_node)
                break

        #res = set()
        #print((frontier))
        #for node in frontier:
        #    res |= node.leaves()

        return frontier

    def query_blacklist(self, mhquery, threshold, scoring_strat, selection_strat, refinement_strat):
        frontier = set([self])

        blacklist = set()
        while True:
            scoreBefore = scoring_strat(frontier, mhquery, threshold)
            candidates = frontier - blacklist
            if len(candidates) < 0.9*len(frontier):
                break
            bad_node = selection_strat(candidates, mhquery, threshold)
            assert(bad_node not in blacklist)

            if bad_node.level == 0:
                frontier.remove(bad_node)
                refinement = set()
            else:
                frontier.remove(bad_node)
                refinement = set(refinement_strat(bad_node, mhquery, threshold)) 
                frontier |= refinement                

            scoreAfter = scoring_strat(frontier, mhquery, threshold)            

            if scoreBefore > scoreAfter:
                frontier -= refinement
                frontier.add(bad_node)
                blacklist.add(bad_node)

        #res = set()
        #for node in frontier:
        #    res |= node.leaves()

        return frontier

    def query_best_match(self, mhquery):
        best_score = -1
        best_node = None
        for n in self.nodes():
            score = mhquery.compare(n.minhash)
            if score > best_score:
                best_node, best_score = n, score
        return [best_node]

    def query_level(self, mhquery, level):
        res = []
        for n in self.nodes(lambda x: x.level == level):
            score1 = mhquery.compare(n.minhash)
            score2 = n.minhash.compare(mhquery)

            if score1 >= 0.05 or score2 >= 0.05:
                res.append(n)
        return res

    def query_gather_mins(self, mhquery, level, expand=False):
        matches = []
        for n in self.nodes(lambda x: x.level == level):
            score = mhquery.compare(n.minhash)
            if score >= 0.05:
                matches.append((score, n))

        matches.sort(reverse=True, key=itemgetter(0))                
        
        if expand:
            subnodes = set()
            for (_, n) in matches:
                subnodes.update(n.children)

            matches = []
            for n in subnodes:
                score = mhquery.compare(n.minhash)
                matches.append((score, n))

        matches.sort(key=itemgetter(0))

        res = []
        query_mins = set(mhquery.get_mins())
        covered = set()
        for (_, n) in matches:
            nodemins = set(n.minhash.get_mins())
            if (nodemins - covered).intersection(query_mins):
                res.append(n)
                covered.update(nodemins)

        return res        

    def leaves(self, visited = None):
        if visited is None:
            visited = set([self])
        if self.level==0:
            #assert(self.level == 0)
            return set([self.id])
        res = set()
        #print("leaves", self.id, len(self.children), len(visited))
        for c in self.children:
            if c not in visited:
                visited.add(c)
                res |= c.leaves(visited)
        return res

    def leaves2(self):
        if len(self.children) == 0:
            assert(self.level == 0)
            return set([self])
        res = set()
        for c in self.children:
            res |= c.leaves2()
        return res

    def shadow(self, visited=None):
        if visited is None:
            visited = set([self])
        if self.level == 0:
            return set([self.vertex])
        res = set()
        for c in self.children:
            if c not in visited:
                visited.add(c)
                res |= c.shadow(visited)
        return res        

    def nodes(self, select=None):
        if select == None:
            select = lambda x: True
        else:
            try:
                0 in select
                collection = select
                select = lambda x: x.id in collection
            except:
                pass

        visited = set()
        frontier = [self]
        while len(frontier) > 0:
            curr = frontier[0]
            frontier = frontier[1:]
            if curr not in visited:
                if select(curr):
                    yield curr
                visited.add(curr)
                frontier += curr.children
 
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
        with open(path.join(projectpath, "{}.catlas.{}.gxt".format(projectname, radius)), 'w') as f:
            writer = Writer(f, ['vertex', 'level'], [])
            for level in self.bfs():
                for v in level:
                    writer.add_vertex(v.id, v.size, [v.vertex, v.level])

            for level in self.bfs():
                for v in level:
                    for c in v.children:
                        writer.add_edge(v.id, c.id)
            writer.done()

        if not self.minhash:
            return

        # Store hashes separately in .mxt
        with open(path.join(projectpath, "{}.catlas.{}.mxt".format(projectname, radius)), 'w') as f:
            for level in self.bfs():
                for v in level:
                    f.write(str(v.id)+',')
                    f.write(' '.join(map(str,v.minhash.get_mins())))
                    f.write('\n')

    @staticmethod 
    def read(catlas_gxt, catlas_mxt, radius):
        graph, node_attr, _ = Graph.from_gxt(open(catlas_gxt, 'r'))
        if catlas_mxt:
            hashes = VertexDict.from_mxt(open(catlas_mxt, 'r'))
        else:
            hashes = {}

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
            vertex = int(node_attr[v]['vertex'])
            cat_nodes[v] = CAtlas(v, vertex, 0, size, [], hashes.get(v))

        max_level = max(levels.keys())
        assert(len(levels[max_level]) == 1)
        for l in range(1, max_level+1):
            for v in levels[l]:
                size = node_attr[v]['size']
                children = filter( lambda x: node_attr[x]['level'] < l, graph.neighbours(v))
                children = list(map( lambda x: cat_nodes[x], children))
                vertex = node_attr[v]['vertex']
                try:
                    vertex = int(vertex)
                except ValueError:
                    pass
                cat_nodes[v] = CAtlas(v, vertex, l, size, children, hashes.get(v))

        root_id = next(iter(levels[max_level]))
        return cat_nodes[root_id]
