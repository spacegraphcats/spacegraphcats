"Test a mildly real subset of twofoo data."
import pytest
import tempfile
import shutil
import os
import csv
import hashlib
import sqlite3

import sourmash
import screed

from spacegraphcats.search.catlas import CAtlas
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


@pytest.mark.pairing
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
def test_dump_contigs():
    global _tempdir

    conf = utils.relative_file("spacegraphcats/conf/twofoo-short.yaml")
    target = "dump_contigs"
    status = run_snakemake(conf, verbose=True, outdir=_tempdir, extra_args=[target])
    assert status == 0

    assert os.path.exists(f"{_tempdir}/twofoo-short_k31/contigs.fa.gz")


@pytest.mark.pairing
@pytest.mark.dependency(depends=["test_build_and_search"])
def test_label_reads():
    global _tempdir

    conf = utils.relative_file("spacegraphcats/conf/twofoo-short.yaml")
    target = "index_reads"
    status = run_snakemake(conf, verbose=True, outdir=_tempdir, extra_args=[target])
    assert status == 0

    assert os.path.exists(f"{_tempdir}/twofoo-short/twofoo-short.reads.bgz")


@pytest.mark.pairing
@pytest.mark.dependency(depends=["test_label_reads"])
def test_index_reads_paired():
    # dig into some of the read pairing stuff in detail
    global _tempdir
    from spacegraphcats.cdbg import index_reads

    test_reads = f"{_tempdir}/twofoo-short/twofoo-short.reads.bgz"
    output_index = os.path.join(_tempdir, "xxx.reads.index")

    retval = index_reads.main(
        [
            f"{_tempdir}/twofoo-short_k31_r1",
            test_reads,
            output_index,
            "-k",
            "31",
            "--expect-paired",
        ]
    )

    assert retval == 0
    assert os.path.exists(output_index)

    def get_ids_at_offset(offset):
        db = sqlite3.connect(output_index)
        c = db.cursor()

        c.execute("SELECT DISTINCT cdbg_id FROM sequences WHERE offset=?", (offset,))
        id_list = [offset for (offset,) in c]
        return id_list

    # output from scripts/print_offsets_names.py
    # SRR606249.7757341 5707860790
    # SRR606249.7757341 5707860911
    assert get_ids_at_offset(5707860790) == get_ids_at_offset(5707860911)
    # SRR606249.7759257 5707861032 - not paired
    # SRR606249.7757341/1 5707861153
    # SRR606249.7757341/2 5707861276
    assert get_ids_at_offset(5707861153) == get_ids_at_offset(5707861276)
    # SRR606249.985492/1 5707861399 - not paired
    # SRR606249.989382/2 5707861521 - not paired
    assert get_ids_at_offset(5707861399) != get_ids_at_offset(5707861521)
    # SRR606249.1044514 5707861643 - not paired


@pytest.mark.pairing
@pytest.mark.dependency(depends=["test_label_reads"])
def test_index_reads_nopairing():
    # dig into some of the read pairing stuff in detail; here, turn off pairs
    global _tempdir
    from spacegraphcats.cdbg import index_reads

    test_reads = f"{_tempdir}/twofoo-short/twofoo-short.reads.bgz"
    output_index = os.path.join(_tempdir, "xxx.reads.index")

    retval = index_reads.main(
        [
            f"{_tempdir}/twofoo-short_k31_r1",
            test_reads,
            output_index,
            "-k",
            "31",
            "--ignore-paired",
        ]
    )

    assert retval == 0
    assert os.path.exists(output_index)

    def get_ids_at_offset(offset):
        db = sqlite3.connect(output_index)
        c = db.cursor()

        c.execute("SELECT DISTINCT cdbg_id FROM sequences WHERE offset=?", (offset,))
        id_list = [offset for (offset,) in c]
        return id_list

    # output from scripts/print_offsets_names.py
    # SRR606249.7757341 5707860790
    # SRR606249.7757341 5707860911
    assert get_ids_at_offset(5707860790) == get_ids_at_offset(5707860911)
    # SRR606249.7759257 5707861032 - not paired
    # SRR606249.7757341/1 5707861153
    # SRR606249.7757341/2 5707861276
    assert get_ids_at_offset(5707861153) != get_ids_at_offset(5707861276)
    # SRR606249.985492/1 5707861399 - belong to different unitigs
    # SRR606249.989382/2 5707861521 - belong to different unitigs!
    assert get_ids_at_offset(5707861399) != get_ids_at_offset(5707861521)
    # SRR606249.1044514 5707861643 - not paired


@pytest.mark.dependency(depends=["test_build_and_search"])
def test_check_contigs_vs_unitigs():
    global _tempdir

    bcalm_sig = "twofoo-short_k31/bcalm.unitigs.fa.sig"
    bcalm_out = sourmash.load_one_signature(os.path.join(_tempdir, bcalm_sig))

    catlas_sig = "twofoo-short_k31_r1/contigs.sig"
    catlas_out = sourmash.load_one_signature(os.path.join(_tempdir, catlas_sig))

    assert bcalm_out.similarity(catlas_out) == 1.0


@pytest.mark.dependency(depends=["test_build_and_search"])
def test_extract_reads_paired():
    global _tempdir

    conf = utils.relative_file("spacegraphcats/conf/twofoo-short.yaml")
    target = "extract_reads"
    status = run_snakemake(conf, verbose=True, outdir=_tempdir, extra_args=[target])
    assert status == 0

    reads_file = os.path.join(
        _tempdir, "twofoo-short_k31_r1_search_oh0/63.short.fa.gz.cdbg_ids.reads.gz"
    )
    read_names = [r.name for r in screed.open(reads_file)]

    num_paired = 0
    num_single = 0

    last_name = None
    for name in read_names:
        is_paired = False
        if last_name:
            if last_name == name:
                is_paired = True
            elif (
                last_name.endswith("/1")
                and name.endswith("/2")
                and last_name[:-2] == name[:-2]
            ):
                is_paired = True

        if is_paired:
            num_paired += 1
        else:
            num_single += 1

        last_name = name

    print(f"reads: {len(read_names)}; singletons: {num_single}; paired: {num_paired}")
    assert len(read_names) == 988
    assert num_single == 616  # this is singletons + first reads
    assert num_paired == 372

    assert num_single + num_paired == len(read_names)


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

    assert m.hexdigest() == "fc9ee74aa29e5d72d2a08c40eee5a0f4", m.hexdigest()

    assert m2.hexdigest() == "92ca814ba49a72022be556aacda59782", m2.hexdigest()


@pytest.mark.dependency(depends=["test_build_and_search"])
def test_check_catlas_vs_contigs():
    global _tempdir

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
