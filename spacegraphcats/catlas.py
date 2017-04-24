"""Data structure for CAtlas."""

import argparse
import cProfile
from spacegraphcats.rdomset import rdomset, domination_graph
from spacegraphcats.graph_io import read_from_gxt


class CAtlas:
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
    def build(graph, radius, domfile):
        """Build a CAtlas at a given radius."""
        # set up values for the base level
        curr_graph = graph
        prev_nodes = None
        level = 0
        idx = 0
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

            # adjust counters
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

            idx += len(nodes)
            level += 1

            # quit if our level is sufficiently small
            if len(domgraph) <= CAtlas.LEVEL_THRESHOLD or \
                    len(domgraph) == len(curr_graph):
                break

            # otherwise prep for the next iteration
            curr_graph = domgraph
            prev_nodes = nodes
        # create a single root over the top level
        root_children = list(nodes.values())
        root_vertex = root_children[0].vertex
        return CAtlas(idx, root_vertex, level, root_children)

    @staticmethod
    def _build_level(graph, radius, level, min_id=0, prev_nodes=None):
        # find the domgraph of the current domgraph
        domset = rdomset(graph, radius)
        domgraph, closest_dominators = domination_graph(graph, domset, radius)

        # closest_dominators indicates the domset vertices that dominate
        # each vertex.
        # v dominating u indicates that u will be a child of v
        # we have the assignment from vertices to dominators, make the
        # reverse
        dominated = {v: list() for v in domset}
        for u, doms in closest_dominators.items():
            for v in doms:
                dominated[v].append(u)

        # create the CAtlas nodes
        nodes = {}
        for idx, v in enumerate(domset):
            # if no previous nodes were supplied, we assume we are on the
            # bottom level and thus the children field is empty
            if prev_nodes is None:
                children = []
            else:
                children = [prev_nodes[u] for u in dominated[v]]
            nodes[v] = CAtlas(min_id+idx, v, level, children)

        return nodes, domgraph, dominated

    def leaves(self, visited=None):
        """Find the descendants of this node with no children."""
        # this function is recursive so we need to keep track of nodes we
        # already visited
        if visited is None:
            visited = set([self])
        # base case is level 0
        if self.level == 0:
            return set([self])
        # otherwise gather the leaves of the children
        res = set()
        for c in self.children:
            if c not in visited:
                visited.add(c)
                res |= c.leaves(visited)
        return res

    def write(self, outfile):
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


def main(args):
    """Build a CAtlas for the provided input graph."""
    r = args.radius
    print("reading graph")
    G = read_from_gxt(args.input, r, False)
    print("reading complete")
    print("building catlas")
    cat = CAtlas.build(G, r, args.domfile)
    print("writing graph")
    cat.write(args.catlas)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input gxt file",
                        type=argparse.FileType('r'))
    parser.add_argument("catlas", help="Output catlas file",
                        type=argparse.FileType('w'))
    parser.add_argument("domfile", help="file to write base domination set",
                        type=argparse.FileType('w'))
    parser.add_argument("radius", help="Catlas radius", type=int)
    args = parser.parse_args()
    cProfile.run("main(args)")
