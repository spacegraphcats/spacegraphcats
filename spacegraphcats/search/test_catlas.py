import os
from .catlas import CAtlas

thisdir = os.path.dirname(__file__)


class Test_LoadCatlas(object):
    def setup(self):
        catlas_prefix = os.path.join(thisdir, 'test-data/catlas.dory_k21_r1')
        catlas_file = os.path.join(catlas_prefix, 'catlas.csv')
        domfile = os.path.join(catlas_prefix, 'first_doms.txt')
        self.catlas = CAtlas(catlas_file, domfile=domfile)
        
    def test_root(self):
        assert self.catlas.root == 731

    def test_levels(self):
        assert len(self.catlas.levels) == 732

    def test_len(self):
        assert len(self.catlas) == 731

    def test_shadow_sizes(self):
        total_shadow_size = 0
        for node_id, level in self.catlas.levels.items():
            if level == 1:
                total_shadow_size += len(self.catlas.layer1_to_cdbg[node_id])

        assert total_shadow_size == 736

        self.catlas.decorate_with_shadow_sizes()
        assert self.catlas.shadow_sizes[self.catlas.root] == total_shadow_size
