import itertools
from spacegraphcats.rdomset import rdomset, domination_graph

class CAtlas:
    LEVEL_THRESHOLD = 1000
    def __init__(self, idx, vertex, level, children):
        """
        Arguments:
            idx:  Integer identifier of the node.  A CAtlas with n nodes will have ids
                0,1,...,n-1
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
    def build(graph, radius):
        """
        Build a CAtlas at a given radius
        """
        # each component will have a separate hierarchy to be merged at the end
        component_top_level = []
        idx = 0 # start index of the component
        for comp in graph.components():
            top_level = CAtlas._build_component(graph, comp, radius, idx)
            component_top_level.append(top_level)
            # adjust the node index so we don't get overlapping indices in different 
            # components.  This assumes that the parent node indices always exceed that 
            # of their children.
            idx = max(node.idx for node in top_level)+1

        # make a root node
        root_children = itertools.chain(*component_top_level)
        # pick a root vertex arbitrarily. 
        root_vertex = root_children[0].vertex
        root_level = max(node.level for node in root_children)
        return CAtlas(idx, root_vertex, root_level, root_children)

    @staticmethod
    def _build_component(graph, comp, radius, min_id):
        """
        Build catlas nodes for one component of the graph
        """
        # find the dominating set
        curr_domset = rdomset(graph, radius, comp)
        # create a dominating graph
        curr_domgraph, closest_dominators = domination_graph(graph, 
            curr_domset, radius, comp)
        # create Catlas nodes for the first domgraph
        idx = min_id
        curr_nodes = {v:CAtlas(idx+i, v, 0, None) for i,v in enumerate(curr_domset)}
        idx += len(curr_domset)

        level = 1
        # keep creating progressively smaller graphs until we hit the level threshold
        while True:
            # find the domgraph of the current domgraph
            next_domset = rdomset(curr_domgraph, 1, curr_domset)
            next_domgraph, closest_dominators = domination_graph(curr_domgraph, next_domset, 1, curr_domset)

            # closest_dominators indicates the domset vertices that dominate each vertex.
            # v dominating u indicates that u will be a child of v
            # we have the assignment from vertices to dominators, make the reverse
            dominated = {v:list() for v in next_domset}
            for u, doms in closest_dominators.items():
                for v in doms:
                    dominated[v].append(u)

            # create the CAtlas nodes for the next level
            next_nodes = {}
            for v in next_domset:
                children = [curr_nodes[u] for u in dominated[v]]
                next_nodes[v] = CAtlas(idx, v, level, children)

            # quit if our level is sufficiently small
            if len(next_domgraph) <= LEVEL_THRESHOLD:
                return list(curr_nodes.values())
            
            # otherwise prep for the next iteration
            curr_domgraph = next_domgraph
            curr_domset = next_domset
            curr_nodes = next_nodes
            level += 1

    def leaves(self, visited=None):
        """
        Find the descendants of this node with no children.
        """
        # this function is recursive so we need to keep track of nodes we already 
        # visited 
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

