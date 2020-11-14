"Tests snakemake execution via click CLI module."
import pytest
import tempfile
import shutil
import os

from spacegraphcats.__main__ import run_snakemake
from . import pytest_utils as utils

# NOTE re dependencies (@pytest.mark.dependency):
# - These basically duplicate the snakemake dependencies.
# - they're there for convenience, because...
# - ...if these are wrong, the tests will still succeed, they just may
#   do some extra work in some tests & take longer.


def setup_module(m):
    global _tempdir
    _tempdir = tempfile.mkdtemp(prefix="sgc_test")


def teardown_module(m):
    global _tempdir
    try:
        shutil.rmtree(_tempdir, ignore_errors=True)
    except OSError:
        pass


@pytest.mark.dependency()
def test_dory_build_cdbg():
    global _tempdir

    dory_conf = utils.relative_file("spacegraphcats/conf/dory-test.yaml")
    target = "dory_k21/bcalm.unitigs.fa"
    status = run_snakemake(
        dory_conf, verbose=True, outdir=_tempdir, extra_args=[target]
    )
    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=["test_dory_build_cdbg"])
def test_dory_build_contigs():
    global _tempdir

    dory_conf = utils.relative_file("spacegraphcats/conf/dory-test.yaml")
    target = "dory_k21/bcalm.unitigs.db"
    status = run_snakemake(
        dory_conf, verbose=True, outdir=_tempdir, extra_args=[target]
    )
    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=["test_dory_build_contigs"])
def test_dory_build_catlas():
    global _tempdir

    dory_conf = utils.relative_file("spacegraphcats/conf/dory-test.yaml")
    target = "dory_k21_r1/catlas.csv"
    status = run_snakemake(
        dory_conf, verbose=True, outdir=_tempdir, extra_args=[target]
    )
    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=["test_dory_build_catlas"])
def test_dory_build_kmer_index():
    global _tempdir

    dory_conf = utils.relative_file("spacegraphcats/conf/dory-test.yaml")
    target = "dory_k21_r1/contigs.mphf"
    status = run_snakemake(
        dory_conf, verbose=True, outdir=_tempdir, extra_args=[target]
    )
    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=["test_dory_build_kmer_index"])
def test_dory_search():
    global _tempdir

    dory_conf = utils.relative_file("spacegraphcats/conf/dory-test.yaml")
    target = "dory_k21_r1_search_oh0/results.csv"
    status = run_snakemake(
        dory_conf, verbose=True, outdir=_tempdir, extra_args=[target]
    )
    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=["test_dory_build_kmer_index"])
def test_dory_build_reads_index():
    global _tempdir

    dory_conf = utils.relative_file("spacegraphcats/conf/dory-test.yaml")
    target = "dory_k21_r1/reads.bgz.index"
    status = run_snakemake(
        dory_conf, verbose=True, outdir=_tempdir, extra_args=[target]
    )
    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=["test_dory_build_reads_index", "test_dory_search"])
def test_dory_extract_reads():
    global _tempdir

    dory_conf = utils.relative_file("spacegraphcats/conf/dory-test.yaml")
    target = "dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.reads.gz"
    status = run_snakemake(
        dory_conf, verbose=True, outdir=_tempdir, extra_args=[target]
    )
    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))


@pytest.mark.dependency(depends=["test_dory_search"])
def test_dory_extract_contigs():
    global _tempdir

    dory_conf = utils.relative_file("spacegraphcats/conf/dory-test.yaml")
    target = "dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.reads.gz"
    status = run_snakemake(
        dory_conf, verbose=True, outdir=_tempdir, extra_args=[target]
    )
    assert status == 0
    assert os.path.exists(os.path.join(_tempdir, target))
