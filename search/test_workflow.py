import os.path
import tempfile
import shutil

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


def test_simple_tr():
    with TempDirectory() as location:
        projpath = os.path.join(location, 'tr')

        # build cDBG
        walk_args = Args()
        walk_args.ksize = 31
        walk_args.memory = 1e8
        walk_args.seqfiles = [relative_filename('data/tr-cross.fa')]
        walk_args.label = False
        walk_args.output = relative_filename(projpath)
        walk_args.loadgraph = None
        walk_args.force = False
        walk_args.no_assemble = False

        walk_dbg.run(walk_args)

        # build catlas
        catlas_args = Args()
        catlas_args.project = relative_filename(projpath)
        catlas_args.radius = 1
        catlas_args.no_checkpoint = True
        catlas_args.level = None

        catlas.main(catlas_args)

        # build minhashes
        args = [projpath]
        args += '-k 31 --scaled=1000'.split(' ')
        make_catlas_minhashes.main(args)
