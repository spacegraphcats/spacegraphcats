import os.path
from . import graph_parser


class CAtlasReader(object):
    """
    Load the given catlas into 'edges', 'vertices', 'roots', 'leaves',
    and 'leaves_to_domnode'.
    """
    def __init__(self, prefix, radius):
        self.prefix = prefix
        self.radius = radius

        self._load()

    def find_level0(self, node_id):
        """Take a node in the catlas and find all level 0 nodes beneath it.

        Return a set of domination graph node IDs.
        """

        x = set()

        beneath = self.edges.get(node_id, [])
        if beneath:
            # recurse!
            for y in beneath:
                x.update(self.find_level0(y))
            return x
        else:
            # only thing there should be nothing beneath are the leaves...
            assert node_id in self.leaves
            # convert leaves into the original domination graph IDs.
            return set([self.leaves_to_domnode[node_id]])

    def _load(self):
        basename = os.path.basename(self.prefix)
        catlas_gxt = '%s.catlas.%d.gxt' % (basename, self.radius)
        catlas_gxt = os.path.join(self.prefix, catlas_gxt)
        self.catlas_gxt = catlas_gxt

        assignment_vxt = '%s.assignment.%d.vxt' % (basename, self.radius)
        assignment_vxt = os.path.join(self.prefix, assignment_vxt)
        self.assignment_vxt = assignment_vxt

        catlas_mxt = '%s.catlas.%d.mxt' % (basename, self.radius)
        catlas_mxt = os.path.join(self.prefix, catlas_mxt)
        self.catlas_mxt = catlas_mxt

        original_graph = '%s.gxt' % (basename)
        original_graph = os.path.join(self.prefix, original_graph)
        self.original_graph = original_graph

        edges = {}
        vertices = {}
        roots = []
        leaves = []
        leaves_to_domnode = {}

        def add_edge(a, b, *extra):
            x = edges.get(a, [])
            x.append(b)
            edges[a] = x

        def add_vertex(node_id, size, names, vals):
            assert names[0] == 'vertex'
            assert names[1] == 'level'
            vertex, level = vals
            if vertex == 'root':
                roots.append(node_id)

            if level == '0':
                leaves.append(node_id)
                leaves_to_domnode[node_id] = int(vertex)

        graph_parser.parse(open(catlas_gxt), add_vertex, add_edge)

        self.edges = edges
        self.vertices = vertices
        self.roots = roots
        self.leaves = leaves
        self.leaves_to_domnode = leaves_to_domnode
    
