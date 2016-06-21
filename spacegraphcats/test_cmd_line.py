#! /usr/bin/env python3
import sys, os.path
sys.path.insert(0, os.path.dirname(__file__))

from . import spg_tst_utils as utils


def test_walk_dbg():
    # run with no args, first
    status, out, err = utils.runscript('walk-dbg.py', [], fail_ok=True)
    assert status == 2


def test_build_cdbg_tr_small():
    tr_small = utils.get_test_data('tr-small.fa')

    with utils.TempDirectory() as tempdir:
        status, out, err = utils.runscript('walk-dbg.py', [tr_small],
                                           in_directory=tempdir)


def test_build_catlas_tr_cross():
    tr_cross = utils.get_test_data('tr-cross.fa')

    with utils.TempDirectory() as tempdir:
        status, out, err = utils.runscript('walk-dbg.py', [tr_cross],
                                           in_directory=tempdir)

        status, out, err = utils.runscript('build-catlas.py',
                                           ['tr-cross', '5'],
                                           in_directory=tempdir)
        print(out)
        assert 'Catlas done' in out
