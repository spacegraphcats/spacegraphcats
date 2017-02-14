import os.path
from . import graph_parser


class CAtlasReader(object):
    """
    Load the given catlas into 'edges', 'vertices', 'roots', 'levels',
    and 'catlas_id_to_domnode'.
    """
    def __init__(self, prefix, radius):
        self.prefix = prefix
        self.radius = radius

        self._load()

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
