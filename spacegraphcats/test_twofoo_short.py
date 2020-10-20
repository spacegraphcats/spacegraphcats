"Test a mildly real subset of twofoo data."
import pytest
import tempfile
import shutil
import os
import csv
import hashlib

import sourmash
import screed

from .search.catlas import CAtlas
from .click import run_snakemake
from .utils import pytest_utils as utils

# NOTE re dependencies (@pytest.mark.dependency):
# - These basically duplicate the snakemake dependencies.
# - they're there for convenience, because...
# - ...if these are wrong, the tests will still succeed, they just may
#   do some extra work in some tests & take longer.


def setup_module(m):
    global _tempdir
    _tempdir = tempfile.mkdtemp(prefix='sgc_test')


def teardown_module(m):
    global _tempdir
    try:
        shutil.rmtree(_tempdir, ignore_errors=True)
    except OSError:
        pass


@pytest.mark.dependency()
def test_build_and_search():
    global _tempdir

    dory_conf = utils.relative_file('spacegraphcats/conf/twofoo-short.yaml')
    target = 'search'
    status = run_snakemake(dory_conf, verbose=True, outdir=_tempdir,
                           extra_args=[target])
    assert status == 0

    output_files = ['twofoo-short/bcalm.twofoo-short.k31.unitigs.fa',
                    'twofoo-short_k31_r1/catlas.csv',
                    'twofoo-short_k31_r1/contigs.fa.gz',
                    'twofoo-short_k31_r1_search_oh0/results.csv']

    for filename in output_files:
        fullpath = os.path.join(_tempdir, filename)
        assert os.path.exists(fullpath), fullpath


@pytest.mark.dependency(depends=['test_build_and_search'])
def test_check_contigs_vs_unitigs():
    global _tempdir

    bcalm_sig = 'twofoo-short/bcalm.twofoo-short.k31.unitigs.fa.sig'
    bcalm_out = sourmash.load_one_signature(os.path.join(_tempdir, bcalm_sig))

    catlas_sig = 'twofoo-short_k31_r1/contigs.fa.gz.sig'
    catlas_out = sourmash.load_one_signature(os.path.join(_tempdir, catlas_sig))

    assert round(bcalm_out.similarity(catlas_out), 2) == .30, bcalm_out.similarity(catlas_out)


@pytest.mark.dependency(depends=['test_build_and_search'])
def test_check_results():
    global _tempdir

    results_csv = os.path.join(_tempdir, 'twofoo-short_k31_r1_search_oh0/results.csv')

    with open(results_csv, 'rt') as fp:
        r = csv.DictReader(fp)

        d = {}
        for row in r:
            q = os.path.basename(row['query'])
            cont = float(row['best_containment'])
            d[q] = cont

    assert len(d) == 3
    assert d['2.short.fa.gz'] == 0.0
    assert round(d['47.short.fa.gz'], 2) == 0.22, round(d['47.short.fa.gz'], 2)
    assert round(d['63.short.fa.gz'], 2) == 0.34, round(d['63.short.fa.gz'], 2)


@pytest.mark.dependency(depends=['test_build_and_search'])
def test_check_md5():
    global _tempdir

    gxt = os.path.join(_tempdir, 'twofoo-short_k31_r1/cdbg.gxt')
    catlas = os.path.join(_tempdir, 'twofoo-short_k31_r1/catlas.csv')

    with open(gxt, 'rb') as fp:
        data = fp.read()
    m = hashlib.md5()
    m.update(data)
    assert m.hexdigest() == '2da4c9bdc4b8fb4be8d239091a6b9b49', m.hexdigest()

    with open(catlas, 'rb') as fp:
        data = fp.read()
    m = hashlib.md5()
    m.update(data)
    assert m.hexdigest() == 'd872c14927854f83a9c2b9d17180f2b5', m.hexdigest()


@pytest.mark.dependency(depends=['test_build_and_search'])
def test_check_catlas_vs_contigs():
    global _tempdir

    catlas_prefix = os.path.join(_tempdir, 'twofoo-short_k31_r1')
    catlas = CAtlas(catlas_prefix)
    print('loaded {} nodes from catlas {}', len(catlas), catlas_prefix)
    print('loaded {} layer 1 catlas nodes', len(catlas.layer1_to_cdbg))

    root = catlas.root
    root_cdbg_nodes = set(catlas.shadow([ root ]))
    print(f'root cdbg nodes: {len(root_cdbg_nodes)}')

    cdbg_id_set = set()
    for record in screed.open(f'{catlas_prefix}/contigs.fa.gz'):
        cdbg_id_set.add(int(record.name))

    print(f'cdbg ID set: {len(cdbg_id_set)}')

    print(f'root - unitigs: {len(root_cdbg_nodes - cdbg_id_set)}')
    print(f'unitigs - root: {len(cdbg_id_set - root_cdbg_nodes)}')

    # all unitigs should be in root shadow
    assert not cdbg_id_set - root_cdbg_nodes
    # all nodes in root shadow should be in unitigs
    assert not root_cdbg_nodes - cdbg_id_set
