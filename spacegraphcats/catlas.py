"""Data structure for CAtlas."""

import argparse
import cProfile
import os
import tempfile
from .rdomset import rdomset, domination_graph
from .graph_io import read_from_gxt, write_to_gxt
from .graph import Graph
from io import TextIOWrapper
from typing import List, Dict, Set


class CAtlas(object):
    """Hierarchical atlas for querying graphs."""

    LEVEL_THRESHOLD = 10

    def __init__(self, idx, vertex, level, children):
        """
        Construct a CAtlas node.

        Arguments:
            idx:  Integer identifier of the node.  A CAtlas with n nodes will
                  have ids 0,1,...,n-1
            vertex:  Name of vertex in the cDBG
            level:  The height of the node in the hierarchy.  The leaves are at
                    level 0, their parents at level 1, etc.
            children:  the CAtlas nodes for which this is a parent
        """
        self.idx = idx
        self.vertex = vertex
        self.children = children
        self.level = level

    @staticmethod
    def build(graph: Graph, radius: int, domfile: TextIOWrapper, 
              checkpoint = None):
        """Build a CAtlas at a given radius."""
        # set up values for the base level if we don't have a checkpoint saved
        if checkpoint is None or not os.path.exists(checkpoint.infile):
            curr_graph = graph
            prev_nodes = None
            level = 0
            idx = 0
        else:
            prev_nodes, curr_graph, idx, level = checkpoint.load()
        # keep creating progressively smaller graphs until we hit the level
        # threshold or steady state
        while True:
            # the base level should have a large radius, others are just 1
            if level == 0:
                r = radius
            else:
                r = 1
            # build the current level
            nodes, domgraph, dominated = CAtlas._build_level(curr_graph, r,
                                                             level, idx,
                                                             prev_nodes)

            print("Catlas level {} complete".format(level))

            # at the bottom level we need to write out the domination
            # assignment
            if level == 0:
                for v, shadow in dominated.items():
                    domstr = str(v)
                    for u in shadow:
                        domstr += " {}".format(u)
                    domstr += "\n"
                    domfile.write(domstr)

            # quit if our level is sufficiently small
            if len(domgraph) <= CAtlas.LEVEL_THRESHOLD or \
                    len(domgraph) == len(curr_graph):
                break

            # prep for the next iteration
            curr_graph = domgraph
            prev_nodes = nodes
            idx += len(nodes)
            level += 1
            
            # write level results to the checkpoint file if applicable
            if checkpoint is not None:
                checkpoint.save(prev_nodes, curr_graph, idx, level)

        # create a single root over the top level
        root_children = list(nodes.values())
        root_vertex = root_children[0].vertex
        return CAtlas(idx, root_vertex, level, root_children)

    @staticmethod
    def _build_level(graph: Graph, radius: int, level: int, min_id: int=0, prev_nodes: List[int]=None):
        # find the domgraph of the current domgraph
        domset = rdomset(graph, radius)
        domgraph, closest_dominators = domination_graph(graph, domset, radius)

        # closest_dominators indicates the domset vertices that dominate
        # each vertex.
        # v dominating u indicates that u will be a child of v
        # we have the assignment from vertices to dominators, make the
        # reverse
        dominated = {v: list() for v in domset}  # type: Dict[int, List[int]]
        for u, doms in closest_dominators.items():
            for v in doms:
                dominated[v].append(u)

        # create the CAtlas nodes
        nodes = {}
        for idx, v in enumerate(domset):
            # if no previous nodes were supplied, we assume we are on the
            # bottom level and thus the children field is empty
            if prev_nodes is None:
                children = []  # type: List[int]
            else:
                children = [prev_nodes[u] for u in dominated[v]]
            nodes[v] = CAtlas(min_id+idx, v, level, children)

        return nodes, domgraph, dominated

    def leaves(self, visited: Set[object]=None) -> Set[object]:
        """Find the descendants of this node with no children."""
        # this function is recursive so we need to keep track of nodes we
        # already visited
        if visited is None:
            visited = set([self])
        # base case is level 0
        if self.level == 0:
            return set([self])
        # otherwise gather the leaves of the children
        res = set()  # type: Set[object]
        for c in self.children:
            if c not in visited:
                visited.add(c)
                res |= c.leaves(visited)
        return res

    def write(self, outfile: TextIOWrapper):
        """Write the connectivity of the CAtlas to file."""
        # doesn't matter how we traverse the graph, so we use DFS for ease of
        # implementation
        stack = [self]
        seen = set()
        while len(stack) > 0:
            # remove from the stack
            curr = stack.pop()
            # write node information
            child_str = " ".join(str(child.idx) for child in curr.children)
            outfile.write("{},{},{},{}\n".format(curr.idx,
                                                 curr.vertex,
                                                 curr.level,
                                                 child_str))
            # all nodes already seen don't get re-added
            seen.add(curr)
            stack.extend(filter(lambda x: x not in seen, curr.children))

    @classmethod
    def read(cls, catlas_file):
        """Load the catlas Directed Acyclic Graph."""
        children = []
        nodes = []

        # load everything from the catlas file
        for line in catlas_file:
            catlas_node, cdbg_node, level, beneath = line.strip().split(',')

            level = int(level)
            catlas_node = int(catlas_node)
            cdbg_node = int(cdbg_node)

            # extend arrays as necessary
            if len(children) <= catlas_node:
                for i in range(catlas_node - len(children) + 1):
                    children.append([])
                    nodes.append(None)

            # parse out the children
            beneath = beneath.strip()
            if beneath:
                beneath = beneath.split(' ')
                children[catlas_node].extend(map(int, beneath))

            # make the new node with empty children
            node = cls(catlas_node, cdbg_node, level, [])
            nodes[catlas_node] = node

        # update the nodes with pointers to their children
        for i, n in enumerate(nodes):
            for child in children[n.idx]:
               n.children.append(nodes[child])

        return nodes[-1] 


class Checkpoint(object):
    """Checkpoints for partial catlas computation."""

    def __init__(self, infile, outfile=None):
        self.infile = infile
        if outfile is None:
            self.outfile = infile
        else:
            self.outfile = outfile

    def load(self):
        """Read cached information from a partial catlas computation."""
        # the temp file contains catlas and graph information.  To use the
        # readers for catlas and graph, we need to temporarily split them into
        # separate files
        tmpf = tempfile.TemporaryFile(mode='r+')
        with open(self.infile, 'r') as f:
            # read until the end of the catlas
            for line in f:
                if line == "###\n":
                    break
                tmpf.write(line)
            # once we are at the graph section, start reading from there
            graph = read_from_gxt(f, radius=1, directed=False, sequential=False)
            # move back to the beginning of the temporary file and read the
            # catlas
            tmpf.seek(0)
            root = CAtlas.read(tmpf)
            tmpf.close()
            print("Root has children {}".format([i.vertex for i in root.children]))
            nodes = {node.vertex: node for node in root.children}
            idx = root.idx
            level = root.level

            if set(graph.nodes) ^ set(nodes.keys()):
                print(graph.nodes)
                print(list(nodes.keys()))
                raise ValueError("graph should have the same nodes as the previous level")

            return nodes, graph, idx, level

    def save(self, nodes, graph, idx, level):
        """Write out a partial computation."""
        with open(self.outfile, 'w') as f:
            # make a dummy root to write the catlas using catlas.write method
            root = CAtlas(idx, 0, level, nodes.values())
            root.write(f)
            f.write("###\n")
            write_to_gxt(f, graph)


def main(args):
    """Build a CAtlas for the provided input graph."""
    r = args.radius
    proj_dir = args.project
    name = os.path.split(proj_dir)[-1]
    graph_file = open(os.path.join(proj_dir, name+".gxt"), 'r')
    catlas_file = open(os.path.join(proj_dir, name+".catlas"), 'w')
    dom_file = open(os.path.join(proj_dir, name+".domfile"), 'w')
    print("reading graph")
    G = read_from_gxt(graph_file, r, False)
    print("reading complete")
    print("building catlas")
    if args.checkpoint is not None:
        cp = Checkpoint(args.checkpoint, args.write_new)
    else:
        cp = None
    cat = CAtlas.build(G, r, dom_file, cp)
    print("writing graph")
    cat.write(catlas_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("project", help="Project directory",
                        type=str)
    parser.add_argument("radius", help="Catlas radius", type=int)
    parser.add_argument("-c", "--checkpoint", nargs='?', type=str,
                        help="Filename of checkpoint")
    parser.add_argument("-w", "--write_new", nargs='?', type=str,
                        help="Filename to write new checkpoints"
                        "If not invoked, old checkpoint will be overwritten")
    args = parser.parse_args()
    main(args)
    #cProfile.run("main(args)")
