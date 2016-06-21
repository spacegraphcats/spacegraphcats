#! /usr/bin/env python3
import sys, os.path
sys.path.insert(0, os.path.dirname(__file__))
import gzip

from . import spg_tst_utils as utils


def test_walk_dbg():
    # run with no args, first
    status, out, err = utils.runscript('walk-dbg.py', [], fail_ok=True)
    assert status == 2


def test_build_cdbg_tr_small():
    return
    tr_small = utils.get_test_data('tr-small.fa')

    with utils.TempDirectory() as tempdir:
        status, out, err = utils.runscript('walk-dbg.py', [tr_small],
                                           in_directory=tempdir)


def test_build_catlas_tr_cross():
    return
    tr_cross = utils.get_test_data('tr-cross.fa')

    with utils.TempDirectory() as tempdir:
        status, out, err = utils.runscript('walk-dbg.py', [tr_cross],
                                           in_directory=tempdir)

        status, out, err = utils.runscript('build-catlas.py',
                                           ['tr-cross', '5'],
                                           in_directory=tempdir)
        print(out)
        assert 'Catlas done' in out


def test_build_catlas_tr_cross_prebuilt():
    tr_cross_gxt = utils.get_test_data('tr-cross.gxt.gz')
    tr_cross_mxt = utils.get_test_data('tr-cross.mxt.gz')

    with utils.TempDirectory() as tempdir:
        catlas_dir = os.path.join(tempdir, 'tr-cross')
        os.mkdir(catlas_dir)

        with open(os.path.join(catlas_dir, 'tr-cross.gxt'), 'wb') as fp:
            with gzip.open(tr_cross_gxt, 'rb') as infp:
                data = infp.read()
                fp.write(data)

        with open(os.path.join(catlas_dir, 'tr-cross.mxt'), 'wb') as fp:
            with gzip.open(tr_cross_mxt, 'rb') as infp:
                data = infp.read()
                fp.write(data)

        status, out, err = utils.runscript('build-catlas.py',
                                           [catlas_dir, '5'])
        print(out)
        assert 'Catlas done' in out
