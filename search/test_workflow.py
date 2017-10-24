import os.path
import tempfile
import shutil
import screed

from spacegraphcats import walk_dbg
from spacegraphcats import catlas
from search import make_catlas_minhashes, frontier_search
import sourmash_lib


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
        walk_args.memory = 1e7
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

        # build query minhash
        max_hash = sourmash_lib.MAX_HASH / 1000.
        query_mh = sourmash_lib.MinHash(ksize=31, n=0, max_hash=max_hash)
        for record in screed.open(relative_filename('data/tr-1.fa')):
            query_mh.add_sequence(record.sequence[:5000])

        print(query_mh.get_mins())
        sig = sourmash_lib.signature.SourmashSignature(query_mh)

        with open('query.sig', 'w') as fp:
            sourmash_lib.signature.save_signatures([ sig ], fp)

        # search!
        args = ['query.sig', projpath, '0.0', '--purgatory',
                '--checkfrontier', '--fullstats']
        frontier_search.main(args)
