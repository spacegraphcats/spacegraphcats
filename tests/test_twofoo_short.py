"Test a mildly real subset of twofoo data."
import pytest
import tempfile
import shutil
import os
import csv
import hashlib

import sourmash
import screed

from spacegraphcats.search.catlas import CAtlas
from spacegraphcats.click import run_snakemake
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
def test_build_and_search():
    global _tempdir

    conf = utils.relative_file("spacegraphcats/conf/twofoo-short.yaml")
    target = "search"
    status = run_snakemake(conf, verbose=True, outdir=_tempdir, extra_args=[target])
    assert status == 0

    output_files = [
        "twofoo-short_k31/bcalm.unitigs.fa",
        "twofoo-short_k31_r1/catlas.csv",
        "twofoo-short_k31_r1/contigs.mphf",
        "twofoo-short_k31_r1_search_oh0/results.csv",
    ]

    for filename in output_files:
        fullpath = os.path.join(_tempdir, filename)
        assert os.path.exists(fullpath), fullpath


@pytest.mark.dependency(depends=["test_build_and_search"])
def test_check_contigs_vs_unitigs():
    global _tempdir

    bcalm_sig = "twofoo-short_k31/bcalm.unitigs.fa.sig"
    bcalm_out = sourmash.load_one_signature(os.path.join(_tempdir, bcalm_sig))

    catlas_sig = "twofoo-short_k31_r1/contigs.sig"
    catlas_out = sourmash.load_one_signature(os.path.join(_tempdir, catlas_sig))

    assert bcalm_out.similarity(catlas_out) == 1.0


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
    assert round(d["47.short.fa.gz"], 2) == 0.52, round(d["47.short.fa.gz"], 2)
    assert round(d["63.short.fa.gz"], 2) == 1.0, round(d["63.short.fa.gz"], 2)


@pytest.mark.dependency(depends=["test_build_and_search"])
def test_check_md5():
    global _tempdir

    gxt = os.path.join(_tempdir, "twofoo-short_k31_r1/cdbg.gxt")
    catlas = os.path.join(_tempdir, "twofoo-short_k31_r1/catlas.csv")

    with open(gxt, "rb") as fp:
        data = fp.read()
    m = hashlib.md5()
    m.update(data)

    with open(catlas, "rb") as fp:
        data = fp.read()
    m2 = hashlib.md5()
    m2.update(data)

    print(m.hexdigest())
    print(m2.hexdigest())

    assert m.hexdigest() == "b14f76a96bf4c5ad2d439009b700c399", m.hexdigest()

    assert m2.hexdigest() == "26982100515262d7b1a380b7b3883ba0", m2.hexdigest()


@pytest.mark.dependency(depends=["test_build_and_search"])
def test_check_catlas_vs_contigs():
    global _tempdir
    import sqlite3

    cdbg_prefix = os.path.join(_tempdir, "twofoo-short_k31")
    catlas_prefix = os.path.join(_tempdir, "twofoo-short_k31_r1")
    catlas = CAtlas(catlas_prefix)
    print("loaded {} nodes from catlas {}", len(catlas), catlas_prefix)
    print("loaded {} layer 1 catlas nodes", len(catlas.layer1_to_cdbg))

    root = catlas.root
    root_cdbg_nodes = set(catlas.shadow([root]))
    print(f"root cdbg nodes: {len(root_cdbg_nodes)}")

    cdbg_id_set = set()
    db = sqlite3.connect(f"{cdbg_prefix}/bcalm.unitigs.db")
    cursor = db.cursor()
    cursor.execute("SELECT id FROM sequences")
    for (record_id,) in cursor:
        cdbg_id_set.add(record_id)

    print(f"cdbg ID set: {len(cdbg_id_set)}")

    print(f"root - unitigs: {len(root_cdbg_nodes - cdbg_id_set)}")
    print(f"unitigs - root: {len(cdbg_id_set - root_cdbg_nodes)}")

    # all unitigs should be in root shadow
    assert not cdbg_id_set - root_cdbg_nodes
    # all nodes in root shadow should be in unitigs
    assert not root_cdbg_nodes - cdbg_id_set
