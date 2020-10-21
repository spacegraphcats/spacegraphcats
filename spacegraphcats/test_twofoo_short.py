"Test a mildly real subset of twofoo data."
import pytest
import tempfile
import shutil
import os
import csv

import sourmash

from .click import run_snakemake
from .utils import pytest_utils as utils

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
def test_build_and_search():
    global _tempdir

    dory_conf = utils.relative_file("spacegraphcats/conf/twofoo-short.yaml")
    target = "search"
    status = run_snakemake(
        dory_conf, verbose=True, outdir=_tempdir, extra_args=[target]
    )
    assert status == 0

    output_files = [
        "twofoo-short/bcalm.twofoo-short.k31.unitigs.fa",
        "twofoo-short_k31_r1/catlas.csv",
        "twofoo-short_k31_r1/contigs.fa.gz",
        "twofoo-short_k31_r1_search_oh0/results.csv",
    ]

    for filename in output_files:
        fullpath = os.path.join(_tempdir, filename)
        assert os.path.exists(fullpath), fullpath


@pytest.mark.dependency(depends=["test_build_and_search"])
def test_check_contigs_vs_unitigs():
    global _tempdir

    bcalm_sig = "twofoo-short/bcalm.twofoo-short.k31.unitigs.fa.sig"
    bcalm_out = sourmash.load_one_signature(os.path.join(_tempdir, bcalm_sig))

    catlas_sig = "twofoo-short_k31_r1/contigs.fa.gz.sig"
    catlas_out = sourmash.load_one_signature(os.path.join(_tempdir, catlas_sig))

    assert round(bcalm_out.similarity(catlas_out), 2) == 0.30, bcalm_out.similarity(
        catlas_out
    )


@pytest.mark.dependency(depends=["test_build_and_search"])
def test_check_results():
    global _tempdir

    results_csv = os.path.join(_tempdir, "twofoo-short_k31_r1_search_oh0/results.csv")

    with open(results_csv, "rt") as fp:
        r = csv.DictReader(fp)

        d = {}
        for row in r:
            q = os.path.basename(row["query"])
            cont = float(row["best_containment"])
            d[q] = cont

    assert len(d) == 3
    assert d["2.short.fa.gz"] == 0.0
    assert round(d["47.short.fa.gz"], 2) == 0.22, round(d["47.short.fa.gz"], 2)
    assert round(d["63.short.fa.gz"], 2) == 0.34, round(d["63.short.fa.gz"], 2)
