import os.path
import tempfile
import shutil
import screed

from spacegraphcats import catlas
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
    newpath = os.path.join(pkgdir, filename)
    return os.path.abspath(newpath)


class Args(object):
    pass


def test_simple_tr():
    with TempDirectory() as location:
        return
    
        projpath = os.path.join(location, 'tr')

        import imp

        # For illustrative purposes.
        runscript = relative_filename('conf/run')
        module = imp.load_source('run', runscript)

        run_args = Args()
        run_args.configfile = 'dory-test'
        run_args.target = 'build'
        run_args.overhead = 0
        run_args.experiment = None
        run_args.radius = 1
        run_args.nolock = False
        run_args.dry_run = False

        module.main(run_args)
