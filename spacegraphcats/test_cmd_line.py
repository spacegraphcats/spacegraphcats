#! /usr/bin/env python3
import sys, os.path
sys.path.insert(0, os.path.dirname(__file__))
import gzip
import shutil

from . import spg_tst_utils as utils


def test_gxt_to_gml():
    # can we run gxt_to_gml?
    tr_cross_gxt = utils.get_test_data('tr-cross.gxt.gz')

    with utils.TempDirectory() as tempdir:
        with open(os.path.join(tempdir, 'tr-cross.gxt'), 'wb') as fp:
            with gzip.open(tr_cross_gxt, 'rb') as infp:
                data = infp.read()
                fp.write(data)

        utils.runscript('gxt_to_gml.py', ['tr-cross.gxt', 'out.gml'],
                        in_directory=tempdir)


def test_walk_dbg():
    # run with no args, first
    status, out, err = utils.runscript('walk-dbg.py', [], fail_ok=True)
    assert status == 2


def test_build_cdbg_tr_small():
    # check a basic run of walk-dbg
    tr_small = utils.get_test_data('tr-small.fa')

    with utils.TempDirectory() as tempdir:
        status, out, err = utils.runscript('walk-dbg.py', [tr_small],
                                           in_directory=tempdir)


def test_build_catlas_tr_cross():
    # can we build a catlas from the output of walk-dbg?
    tr_cross = utils.get_test_data('tr-cross.fa')

    with utils.TempDirectory() as tempdir:
        status, out, err = utils.runscript('walk-dbg.py', [tr_cross],
                                           in_directory=tempdir)

        status, out, err = utils.runscript('build-catlas.py',
                                           ['tr-cross', '5'],
                                           in_directory=tempdir)
        print(out)
        assert 'Catlas done' in out


def _build_catlas_from_gz(gxtfile, mxtfile, tempdir, level=5):
    # construct catlas name from first component of gxtfile
    catlas_name = os.path.basename(gxtfile).split('.')[0]
    catlas_dir = os.path.join(tempdir, catlas_name)
    os.mkdir(catlas_dir)

    # copy over the gz gxt file
    with open(os.path.join(catlas_dir, catlas_name + '.gxt'), 'wb') as fp:
        with gzip.open(gxtfile, 'rb') as infp:
            data = infp.read()
            fp.write(data)

    # copy over the gz mxt file
    with open(os.path.join(catlas_dir, catlas_name + '.mxt'), 'wb') as fp:
        with gzip.open(mxtfile, 'rb') as infp:
            data = infp.read()
            fp.write(data)

    # run build-catlas!
    status, out, err = utils.runscript('build-catlas.py',
                                       [catlas_dir, str(level)])
    print(out)
    # make sure it worked...
    assert 'Catlas done' in out

    return catlas_dir


def test_build_catlas_tr_cross_prebuilt():
    # can we build a catlas from prebuilt gxt/mxt?
    tr_cross_gxt = utils.get_test_data('tr-cross.gxt.gz')
    tr_cross_mxt = utils.get_test_data('tr-cross.mxt.gz')

    with utils.TempDirectory() as tempdir:
        dirname = _build_catlas_from_gz(tr_cross_gxt, tr_cross_mxt,
                                        tempdir)


def test_search_new_catlas():
    # build a catlas and then search it --
    tr_cross_gxt = utils.get_test_data('tr-cross.gxt.gz')
    tr_cross_mxt = utils.get_test_data('tr-cross.mxt.gz')

    mh1_txt = utils.get_test_data('tr-1.fa.sig.dump.txt')
    mh2_txt = utils.get_test_data('tr-2.fa.sig.dump.txt')
    mh3_txt = utils.get_test_data('acido-short.fa.sig.dump.txt')

    with utils.TempDirectory() as tempdir:
        dirname = _build_catlas_from_gz(tr_cross_gxt, tr_cross_mxt,
                                        tempdir)

        status, out, err = utils.runscript('search-catlas-mh.py',
                            [dirname, mh1_txt], in_directory=tempdir)
        assert 'tr-1.fa.sig.dump.txt,0.488,0.976,0,1' in out

        status, out, err = utils.runscript('search-catlas-mh.py',
                            [dirname, mh2_txt], in_directory=tempdir)
        assert 'tr-2.fa.sig.dump.txt,0.500,1.000,0,1' in out

        status, out, err = utils.runscript('search-catlas-mh.py',
                            [dirname, mh3_txt], in_directory=tempdir)

        assert out.endswith('levels\n---\n'), (out,)
        print(out)


def test_search_prebuilt_catlas():
    # use a pre-built catlas and search it --
    tr_cross_dir = utils.get_test_data('catlas.tr-cross')

    mh1_txt = utils.get_test_data('tr-1.fa.sig.dump.txt')
    mh2_txt = utils.get_test_data('tr-2.fa.sig.dump.txt')
    mh3_txt = utils.get_test_data('acido-short.fa.sig.dump.txt')

    with utils.TempDirectory() as tempdir:
        dirname = os.path.join(tempdir, 'tr-cross')
        shutil.copytree(tr_cross_dir, dirname)

        # results
        status, out, err = utils.runscript('search-catlas-mh.py',
                            [dirname, mh1_txt], in_directory=tempdir)
        assert 'tr-1.fa.sig.dump.txt,0.488,0.976,0,1' in out

        # results
        status, out, err = utils.runscript('search-catlas-mh.py',
                            [dirname, mh2_txt], in_directory=tempdir)
        assert 'tr-2.fa.sig.dump.txt,0.500,1.000,0,1' in out

        # no results
        status, out, err = utils.runscript('search-catlas-mh.py',
                            [dirname, mh3_txt], in_directory=tempdir)

        assert out.endswith('levels\n---\n'), (out,)
        print(out)
