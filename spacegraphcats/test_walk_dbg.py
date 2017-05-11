import os.path
import tempfile
import shutil
import screed

from spacegraphcats import walk_dbg
from spacegraphcats import catlas
from search import make_catlas_minhashes


class TempDirectory(object):
    def __init__(self):
        self.tempdir = tempfile.mkdtemp(prefix='sgc_test')

    def __enter__(self):
        return self.tempdir

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            shutil.rmtree(self.tempdir, ignore_errors=True)
        except OSError:
            pass

        if exc_type:
            return False


def relative_filename(filename):
    thisdir = os.path.dirname(__file__)
    pkgdir = os.path.join(thisdir, '..')
    return os.path.join(pkgdir, filename)


class Args(object):
    pass


# this test looks at a purely linear segment in tr_vsmall, and ensures that
# there is a single node and single contig output; since the primary behavior
# of walk_dbg is to traverse from all high-degree nodes, earlier versions
# omitted purely linear paths. Never again!

def test_purely_linear_tr_vsmall():
    with TempDirectory() as location:
        projpath = os.path.join(location, 'tr_vsmall')

        # build cDBG
        walk_args = Args()
        walk_args.ksize = 31
        walk_args.memory = 1e8
        walk_args.seqfiles = [relative_filename('data/tr-vsmall.fa')]
        walk_args.label = False
        walk_args.output = relative_filename(projpath)
        walk_args.loadgraph = None
        walk_args.force = False
        walk_args.no_assemble = False

        walk_dbg.run(walk_args)

        # ok, now look at output --

        cdbg_file = os.path.join(projpath, 'cdbg.gxt')
        num_nodes = len(list(open(cdbg_file)))
        assert num_nodes == 1
        
        contigs_file = os.path.join(projpath, 'contigs.fa.gz')
        num_contigs = len(list(screed.open(contigs_file)))
        assert num_contigs == 1
