"Utilities for pytest-based tests; see test_* under spacegraphcats package."
import shutil
import os
import tempfile


class TempDirectory(object):
    "Set up a temp directory; optionally change into it. Context manager."

    def __init__(self, chdir=False):
        self.tempdir = tempfile.mkdtemp(prefix="sgc_test")
        self.start_directory = None
        self.do_chdir = chdir

    def __enter__(self):
        if self.do_chdir:
            self.start_directory = os.getcwd()
            os.chdir(self.tempdir)
        return self.tempdir

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            shutil.rmtree(self.tempdir, ignore_errors=True)
        except OSError:
            pass

        if self.start_directory:
            os.chdir(self.start_directory)

        if exc_type:
            return False


def relative_file(filename):
    "Return the filename relative to the top level directory of this repo."
    thisdir = os.path.dirname(__file__)
    pkgdir = os.path.join(thisdir, "..")
    newpath = os.path.join(pkgdir, filename)
    return os.path.abspath(newpath)


def pkg_file(filename):
    "Return the filename relative to the spacegraphcats/ package."
    thisdir = os.path.dirname(__file__)
    pkgdir = os.path.join(thisdir, "..")
    newpath = os.path.join(pkgdir, filename)
    return os.path.abspath(newpath)


def in_tempdir(fn):
    "Decorator: run this test function in tempdir, provide location."

    def wrapper(*args, **kwargs):
        with TempDirectory(chdir=True) as location:
            newargs = [location] + list(args)
            return fn(*newargs, **kwargs)

    return wrapper


def in_thisdir(fn):
    "Decorator: run this test function in cwd, provide location of temp dir"

    def wrapper(*args, **kwargs):
        with TempDirectory() as location:
            newargs = [location] + list(args)
            return fn(*newargs, **kwargs)

    return wrapper


class Args(object):
    "Empty object to use for stashing argparse-like args."
    pass
