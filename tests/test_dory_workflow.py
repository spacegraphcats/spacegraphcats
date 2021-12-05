"""
Test many (most?) spacegraphcats scripts against the 'dory' small example data.
"""
import os.path
import shutil
import glob
import hashlib
import gzip
import sqlite3

import screed
import sourmash

from . import pytest_utils
from .pytest_utils import pkg_file, relative_file

from spacegraphcats.catlas import catlas
from spacegraphcats.cdbg import index_cdbg_by_kmer, MPHF_KmerIndex
from spacegraphcats.search import query_by_sequence
from spacegraphcats.search import characterize_catlas_regions
from spacegraphcats.search import extract_unassembled_nodes
from spacegraphcats.search import evaluate_overhead
from spacegraphcats.search import catlas_info
from spacegraphcats.search import extract_contigs
from spacegraphcats.search import extract_contigs_cdbg
from spacegraphcats.search import estimate_query_abundance
from spacegraphcats.search import extract_nodes_by_shadow_ratio
from spacegraphcats.utils import make_bgzf
from spacegraphcats.cdbg import index_reads
from spacegraphcats.search import extract_reads
from spacegraphcats.cdbg import index_cdbg_by_minhash
from spacegraphcats.search import query_by_hashval
from spacegraphcats.search import index_cdbg_by_multifasta
from spacegraphcats.search import index_cdbg_by_multifasta_x
from spacegraphcats.search import query_multifasta_by_sig
from spacegraphcats.search import extract_cdbg_by_multifasta
from spacegraphcats.search import search_utils
from spacegraphcats.search import query_by_prot
from spacegraphcats.search import extract_neighborhoods_by_cdbg_ids
from spacegraphcats.search import count_dominator_abundance


def copy_dory_catlas():
    testdata = pkg_file("tests/test-data/catlas.dory_k21_r1")
    shutil.copytree(testdata, "./dory_k21_r1")

    testdata = pkg_file("tests/test-data/catlas.dory_k21")
    shutil.copytree(testdata, "./dory_k21")


def copy_dory_catlas_search():
    testdata = pkg_file("tests/test-data/catlas.dory_k21_r1_search_oh0")
    shutil.copytree(testdata, "./dory_k21_r1_search_oh0")


def copy_dory_head():
    testdata = relative_file("data/dory-head.fa")
    shutil.copyfile(testdata, "dory-head.fa")


def copy_dory_subset():
    testdata = relative_file("data/dory-subset.fa")
    shutil.copyfile(testdata, "dory-subset.fa")

    testdata = relative_file("data/dory-subset.fq")
    shutil.copyfile(testdata, "dory-subset.fq")


def copy_dory_sig():
    testdata = relative_file("data/dory-subset.fq.sig")
    shutil.copyfile(testdata, "dory-subset.fq.sig")


### actual tests


@pytest_utils.in_tempdir
def test_dory_break_bad_bcalm(location):
    from spacegraphcats.cdbg import sort_bcalm_unitigs

    copy_dory_head()
    copy_dory_subset()

    # make the output directory
    try:
        os.mkdir("dory_k21")
    except FileExistsError:
        pass

    # sort the bcalm file
    args = [
        "-k",
        "21",
        relative_file("data/bcalm-BROKEN.dory.k21.unitigs.fa"),
        "dory_k21/bcalm.unitigs.db",
        "dory_k21/bcalm.unitigs.pickle",
    ]

    assert sort_bcalm_unitigs.main(args) != 0


@pytest_utils.in_tempdir
def test_dory_sort_bcalm_diff_seed(location):
    from spacegraphcats.cdbg import sort_bcalm_unitigs

    copy_dory_head()
    copy_dory_subset()

    # make the output directory
    try:
        os.mkdir("dory_k21_seedtest")
    except FileExistsError:
        pass

    # sort the bcalm file
    args = [
        "-k",
        "21",
        relative_file("data/bcalm-BROKEN.dory.k21.unitigs.fa"),
        "dory_k21_seedtest/bcalm.unitigs.db",
        "dory_k21_seedtest/bcalm.unitigs.pickle",
        "--seed", "43",
    ]

    assert sort_bcalm_unitigs.main(args) != 0


@pytest_utils.in_tempdir
def test_dory_query_workflow(location):
    from spacegraphcats.cdbg import bcalm_to_gxt, sort_bcalm_unitigs

    copy_dory_head()
    copy_dory_subset()

    # make the output directory
    try:
        os.mkdir("dory_k21")
        os.mkdir("dory_k21_r1")
    except FileExistsError:
        pass

    # sort the bcalm file
    args = [
        "-k",
        "21",
        relative_file("data/bcalm.dory.k21.unitigs.fa"),
        "dory_k21/bcalm.unitigs.db",
        "dory_k21/bcalm.unitigs.pickle",
    ]

    assert sort_bcalm_unitigs.main(args) == 0

    db = sqlite3.connect("dory_k21/bcalm.unitigs.db")
    all_seqs = list(search_utils.contigs_iter_sqlite(db))
    assert len(all_seqs) == 736, len(all_seqs)

    # convert the bcalm file to gxt
    args = [
        "-P",
        "dory_k21/bcalm.unitigs.db",
        "dory_k21/bcalm.unitigs.pickle",
        "dory_k21/cdbg.gxt",
        "dory_k21/contigs",
    ]

    assert bcalm_to_gxt.main(args) == 0

    db = sqlite3.connect("dory_k21/bcalm.unitigs.db")
    all_seqs = list(search_utils.contigs_iter_sqlite(db))
    assert len(all_seqs) == 736, len(all_seqs)

    with open("dory_k21/cdbg.gxt", "rb") as fp:
        data = fp.read()
    m = hashlib.md5()
    m.update(data)
    assert m.hexdigest() == "79d74806263900d54078eced695ec193", m.hexdigest()

    # build catlas
    args = pytest_utils.Args()
    args.no_checkpoint = True
    args.level = 0
    args.radius = 1
    args.cdbg_dir = "dory_k21"
    args.catlas_dir = "dory_k21_r1"
    print("** running catlas")
    assert catlas.main(args) == 0

    # make k-mer search index
    args = "-k 21 dory_k21 --contigs-db dory_k21/bcalm.unitigs.db".split()
    print("** running index_cdbg_by_kmer")
    assert index_cdbg_by_kmer.main(args) == 0

    # check that we get the kmer -> cDBG assignments we expect
    kmer_idx = MPHF_KmerIndex.from_directory("dory_k21")
    assert kmer_idx.table[10218271035842461694] == 118
    assert kmer_idx.table[8436068710919520258] == 118
    assert kmer_idx.table[13994045974119358468] == 118
    assert kmer_idx.table[11971930231572094512] == 187

    # do search!!
    args = "dory_k21 dory_k21_r1 dory_k21_r1_search_oh0 --query dory-head.fa -k 21 --contigs-db dory_k21/bcalm.unitigs.db".split()
    try:
        assert query_by_sequence.main(args) == 0
    except SystemExit as e:
        assert e.code == 0, str(e)

    # check output!
    output_path = "dory_k21_r1_search_oh0/"
    assert os.path.exists(output_path + "command.txt")
    assert os.path.exists(output_path + "dory-head.fa.frontier.txt.gz")
    assert os.path.exists(output_path + "dory-head.fa.cdbg_ids.txt.gz")
    assert os.path.exists(output_path + "dory-head.fa.response.txt")
    assert os.path.exists(output_path + "dory-head.fa.contigs.sig")
    assert os.path.exists(output_path + "results.csv")

    with gzip.open(output_path + "dory-head.fa.cdbg_ids.txt.gz", "rt") as fp:
        match_ids = [x.strip() for x in fp]
        assert match_ids == ["118", "187"]

    with open(output_path + "results.csv", "rt") as fp:
        lines = fp.readlines()
        assert len(lines) == 2

        last_line = lines[-1].strip()
        assert (
            last_line == "dory-head.fa,1.0,1.0,1671,2,21,1631,1.0,0.0,0.0,dory_k21_r1"
        )

    sigfile = output_path + "dory-head.fa.contigs.sig"
    with open(sigfile, "rt") as fp:
        sig = sourmash.load_one_signature(fp)
        name = str(sig)
        assert "from dory_k21_r1" in name
        assert "nbhd:TRINITY_DN290219_c0_g1_i1" in name


# CTB note: this test can be merged with test_dory_query_workflow if we
# create checkpoint fixture that supports setting checkpoint to both
# true and false. We'll also need to make the tempdir stuff a fixture, too;
# sourmash has both of these so we can swipe code from there.
@pytest_utils.in_tempdir
def test_dory_query_workflow_checkpoint(location):
    from spacegraphcats.cdbg import bcalm_to_gxt, sort_bcalm_unitigs

    copy_dory_head()
    copy_dory_subset()

    # make the output directory
    try:
        os.mkdir("dory_k21")
        os.mkdir("dory_k21_r1")
    except FileExistsError:
        pass

    # sort the bcalm file
    args = [
        "-k",
        "21",
        relative_file("data/bcalm.dory.k21.unitigs.fa"),
        "dory_k21/bcalm.unitigs.db",
        "dory_k21/bcalm.unitigs.pickle",
    ]

    assert sort_bcalm_unitigs.main(args) == 0

    db = sqlite3.connect("dory_k21/bcalm.unitigs.db")
    all_seqs = list(search_utils.contigs_iter_sqlite(db))
    assert len(all_seqs) == 736, len(all_seqs)

    # convert the bcalm file to gxt
    args = [
        "-P",
        "dory_k21/bcalm.unitigs.db",
        "dory_k21/bcalm.unitigs.pickle",
        "dory_k21/cdbg.gxt",
        "dory_k21/contigs",
    ]

    assert bcalm_to_gxt.main(args) == 0

    db = sqlite3.connect("dory_k21/bcalm.unitigs.db")
    all_seqs = list(search_utils.contigs_iter_sqlite(db))
    assert len(all_seqs) == 736, len(all_seqs)

    with open("dory_k21/cdbg.gxt", "rb") as fp:
        data = fp.read()
    m = hashlib.md5()
    m.update(data)
    assert m.hexdigest() == "79d74806263900d54078eced695ec193", m.hexdigest()

    # build catlas
    args = pytest_utils.Args()
    args.no_checkpoint = False
    args.level = 0
    args.radius = 1
    args.cdbg_dir = "dory_k21"
    args.catlas_dir = "dory_k21_r1"
    print("** running catlas")
    assert catlas.main(args) == 0

    # make k-mer search index
    args = "-k 21 dory_k21 --contigs-db dory_k21/bcalm.unitigs.db".split()
    print("** running index_cdbg_by_kmer")
    assert index_cdbg_by_kmer.main(args) == 0

    # check that we get the kmer -> cDBG assignments we expect
    kmer_idx = MPHF_KmerIndex.from_directory("dory_k21")
    assert kmer_idx.table[10218271035842461694] == 118
    assert kmer_idx.table[8436068710919520258] == 118
    assert kmer_idx.table[13994045974119358468] == 118
    assert kmer_idx.table[11971930231572094512] == 187

    # do search!!
    args = "dory_k21 dory_k21_r1 dory_k21_r1_search_oh0 --query dory-head.fa -k 21 --contigs-db dory_k21/bcalm.unitigs.db".split()
    try:
        assert query_by_sequence.main(args) == 0
    except SystemExit as e:
        assert e.code == 0, str(e)

    # check output!
    output_path = "dory_k21_r1_search_oh0/"
    assert os.path.exists(output_path + "command.txt")
    assert os.path.exists(output_path + "dory-head.fa.frontier.txt.gz")
    assert os.path.exists(output_path + "dory-head.fa.cdbg_ids.txt.gz")
    assert os.path.exists(output_path + "dory-head.fa.response.txt")
    assert os.path.exists(output_path + "dory-head.fa.contigs.sig")
    assert os.path.exists(output_path + "results.csv")

    with gzip.open(output_path + "dory-head.fa.cdbg_ids.txt.gz", "rt") as fp:
        match_ids = [x.strip() for x in fp]
        assert match_ids == ["118", "187"]

    with open(output_path + "results.csv", "rt") as fp:
        lines = fp.readlines()
        assert len(lines) == 2

        last_line = lines[-1].strip()
        assert (
            last_line == "dory-head.fa,1.0,1.0,1671,2,21,1631,1.0,0.0,0.0,dory_k21_r1"
        )

    sigfile = output_path + "dory-head.fa.contigs.sig"
    with open(sigfile, "rt") as fp:
        sig = sourmash.load_one_signature(fp)
        name = str(sig)
        assert "from dory_k21_r1" in name
        assert "nbhd:TRINITY_DN290219_c0_g1_i1" in name


@pytest_utils.in_tempdir
def test_dory_query_workflow_remove_pendants(location):
    from spacegraphcats.cdbg import bcalm_to_gxt, sort_bcalm_unitigs

    copy_dory_head()
    copy_dory_subset()

    # make the output directory
    try:
        os.mkdir("dory_k21")
        os.mkdir("dory_k21_r1")
    except FileExistsError:
        pass

    # sort the bcalm file
    args = [
        "-k",
        "21",
        relative_file("data/bcalm.dory.k21.unitigs.fa"),
        "dory_k21/bcalm.unitigs.db",
        "dory_k21/bcalm.unitigs.pickle",
    ]

    assert sort_bcalm_unitigs.main(args) == 0

    db = sqlite3.connect("dory_k21/bcalm.unitigs.db")
    all_seqs = list(search_utils.contigs_iter_sqlite(db))
    assert len(all_seqs) == 736, len(all_seqs)

    # convert the bcalm file to gxt
    args = [
        "dory_k21/bcalm.unitigs.db",
        "dory_k21/bcalm.unitigs.pickle",
        "dory_k21/cdbg.gxt",
        "dory_k21/contigs",
    ]

    assert bcalm_to_gxt.main(args) == 0

    db = sqlite3.connect("dory_k21/bcalm.unitigs.db")
    all_seqs = list(search_utils.contigs_iter_sqlite(db))
    assert len(all_seqs) == 736, len(all_seqs)

    with open("dory_k21/cdbg.gxt", "rb") as fp:
        data = fp.read()
    m = hashlib.md5()
    m.update(data)
    assert m.hexdigest() == "7e4d9acc9e968f7425c94f6ec78ecdd5", m.hexdigest()


@pytest_utils.in_tempdir
def test_dory_search_nomatch(location):
    # test situations where zero k-mers match - should not fail.
    copy_dory_catlas()

    testdata = relative_file("data/random-query-nomatch.fa")
    shutil.copyfile(testdata, "random-query.fa")

    # make k-mer search index
    args = "-k 21 dory_k21 --contigs-db dory_k21/bcalm.unitigs.db".split()
    print("** running index_cdbg_by_kmer")
    assert index_cdbg_by_kmer.main(args) == 0

    # do search!!
    args = "dory_k21 dory_k21_r1 dory_k21_r1_search_oh0 --query random-query.fa -k 21 --contigs-db dory_k21/bcalm.unitigs.db".split()
    try:
        assert query_by_sequence.main(args) == 0
    except SystemExit as e:
        assert e.code == 0, str(e)


@pytest_utils.in_tempdir
def test_dory_characterize_catlas_regions(location):
    copy_dory_catlas()
    copy_dory_head()
    copy_dory_subset()

    # run characterize_catlas_regions
    args = "dory_k21 dory_k21_r1 dory_k1_r1.vec --contigs-db dory_k21/bcalm.unitigs.db".split()
    assert characterize_catlas_regions.main(args) == 0


@pytest_utils.in_tempdir
def test_dory_extract_unassembled_nodes(location):
    copy_dory_catlas()
    copy_dory_head()

    # make k-mer search index - FIXTURE
    args = "-k 21 dory_k21 --contigs-db dory_k21/bcalm.unitigs.db".split()
    assert index_cdbg_by_kmer.main(args) == 0

    # run extract_unassembled_regions
    args = "dory_k21 dory_k21_r1 dory-head.fa dory.regions -k 21 --contigs-db dory_k21/bcalm.unitigs.db".split()
    assert extract_unassembled_nodes.main(args) == 0


@pytest_utils.in_tempdir
def test_dory_catlas_info(location):
    copy_dory_catlas()

    # run catlas info
    print("running catlas_info")
    assert catlas_info.main(["dory_k21", "dory_k21_r1"]) == 0


@pytest_utils.in_tempdir
def test_dory_extract_contigs(location):
    copy_dory_catlas()
    copy_dory_catlas_search()
    copy_dory_head()

    # run extract_contigs
    print("running extract_contigs")
    args = [
        "--contigs-db",
        "dory_k21/bcalm.unitigs.db",
        "dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.txt.gz",
        "-o",
        "dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.contigs.fa.gz",
    ]
    assert extract_contigs.main(args) == 0

    assert os.path.exists("dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.contigs.fa.gz")


@pytest_utils.in_tempdir
def test_dory_extract_contigs_cdbg(location):
    copy_dory_catlas()
    copy_dory_catlas_search()
    copy_dory_head()

    # make k-mer search index - FIXTURE
    args = "-k 21 dory_k21_r1 --contigs-db dory_k21/bcalm.unitigs.db".split()
    assert index_cdbg_by_kmer.main(args) == 0

    # run extract_contigs
    print("running extract_contigs_cdbg")
    args = [
        "dory_k21_r1",
        "dory-head.fa",
        "--contigs-db",
        "dory_k21/bcalm.unitigs.db",
        "-o",
        "dory_k21_r1_search_oh0/dory-head.fa.matches.contigs.fa.gz",
        "-k",
        "21",
    ]
    assert extract_contigs_cdbg.main(args) == 0

    assert os.path.exists("dory_k21_r1_search_oh0/dory-head.fa.matches.contigs.fa.gz")


@pytest_utils.in_tempdir
def test_dory_make_bgzf(location):
    copy_dory_subset()

    # run make_bgzf
    print("** running make_bgzf")
    args = ["dory-subset.fa", "-o", "reads.bgz"]
    assert make_bgzf.main(args) == 0


@pytest_utils.in_tempdir
def test_dory_index_reads(location):
    copy_dory_catlas()
    copy_dory_subset()

    # run make_bgzf - FIXTURE
    print("** running make_bgzf")
    args = ["dory-subset.fa", "-o", "reads.bgz"]
    assert make_bgzf.main(args) == 0

    # make k-mer search index - FIXTURE
    args = "-k 21 dory_k21_r1 --contigs-db dory_k21/bcalm.unitigs.db".split()
    assert index_cdbg_by_kmer.main(args) == 0

    # run index_reads
    print("** running index_reads")
    args = [
        "-k",
        "21",
        "dory_k21_r1",
        "reads.bgz",
        "dory_k21_r1/reads.bgz.index",
    ]
    assert index_reads.main(args) == 0


@pytest_utils.in_tempdir
def test_dory_index_reads_require_paired_fail(location):
    copy_dory_catlas()
    copy_dory_subset()

    # run make_bgzf - FIXTURE
    print("** running make_bgzf")
    args = ["dory-subset.fa", "-o", "reads.bgz"]
    assert make_bgzf.main(args) == 0

    # make k-mer search index - FIXTURE
    args = "-k 21 dory_k21_r1 --contigs-db dory_k21/bcalm.unitigs.db".split()
    assert index_cdbg_by_kmer.main(args) == 0

    # run index_reads
    print("** running index_reads")
    args = [
        "-k",
        "21",
        "dory_k21_r1",
        "reads.bgz",
        "dory_k21_r1/reads.bgz.index",
        "-P",
    ]
    assert index_reads.main(args) != 0


@pytest_utils.in_tempdir
def test_dory_index_reads_check_args_fail(location):
    # run index_reads
    print("** running index_reads")
    args = [
        "-k",
        "21",
        "dory_k21_r1",
        "reads.bgz",
        "dory_k21_r1/reads.bgz.index",
        "-P",
        "-N",
    ]
    assert index_reads.main(args) != 0


@pytest_utils.in_tempdir
def test_dory_extract_reads(location):
    copy_dory_catlas()
    copy_dory_catlas_search()
    copy_dory_subset()

    # run make_bgzf - FIXTURE
    print("** running make_bgzf")
    args = ["dory-subset.fa", "-o", "reads.bgz"]
    assert make_bgzf.main(args) == 0

    # make k-mer search index - FIXTURE
    args = "-k 21 dory_k21_r1 --contigs-db dory_k21/bcalm.unitigs.db".split()
    assert index_cdbg_by_kmer.main(args) == 0

    # run label_cdbg - FIXTURE
    print("** running label_cdbg")
    args = [
        "-k",
        "21",
        "dory_k21_r1",
        "reads.bgz",
        "dory_k21_r1/reads.bgz.index",
    ]
    assert index_reads.main(args) == 0

    # run extract_reads
    print("** running extract_reads")
    args = [
        "reads.bgz",
        "dory_k21_r1/reads.bgz.index",
        "dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.txt.gz",
        "-o",
        "dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.reads.gz",
    ]
    assert extract_reads.main(args) == 0

    reads_filename = "dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.reads.gz"
    reads = [record for record in screed.open(reads_filename)]
    assert len(reads) == 2


@pytest_utils.in_tempdir
def test_dory_extract_reads_fq(location):
    # check that FASTQ is passed through properly.

    copy_dory_catlas()
    copy_dory_catlas_search()
    copy_dory_subset()

    # run make_bgzf - FIXTURE
    print("** running make_bgzf")
    args = ["dory-subset.fq", "-o", "reads.bgz"]
    assert make_bgzf.main(args) == 0

    # make k-mer search index - FIXTURE
    args = "-k 21 dory_k21_r1 --contigs-db dory_k21/bcalm.unitigs.db".split()
    assert index_cdbg_by_kmer.main(args) == 0

    # run index_reads - FIXTURE
    print("** running label_cdbg")
    args = [
        "-k",
        "21",
        "dory_k21_r1",
        "reads.bgz",
        "dory_k21_r1/reads.bgz.labels",
    ]
    assert index_reads.main(args) == 0

    # run extract_reads
    print("** running extract_reads")
    args = [
        "reads.bgz",
        "dory_k21_r1/reads.bgz.labels",
        "dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.txt.gz",
        "-o",
        "dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.reads.gz",
    ]
    assert extract_reads.main(args) == 0

    reads_filename = "dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.reads.gz"
    reads = [record for record in screed.open(reads_filename)]
    assert len(reads) == 2
    assert len(reads[0].quality)  # FASTQ preserved!


@pytest_utils.in_tempdir
def test_dory_evaluate_overhead(location):
    copy_dory_catlas()
    copy_dory_head()
    copy_dory_catlas_search()

    # run evaluate_overhead
    args = [
        "dory_k21",
        "dory_k21_r1",
        "dory-head.fa",
        "dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.txt.gz",
        "-o",
        "xyz",
        "--contigs-db",
        "dory_k21/bcalm.unitigs.db",
    ]
    print("** running evaluate_overhead")
    assert evaluate_overhead.main(args) == 0


@pytest_utils.in_tempdir
def test_dory_estimate_query_abundance(location):
    copy_dory_catlas()
    copy_dory_head()

    # make k-mer search index - FIXTURE
    args = "-k 21 dory_k21 --contigs-db dory_k21/bcalm.unitigs.db".split()
    assert index_cdbg_by_kmer.main(args) == 0

    # calculate query abundances
    args = "dory_k21 dory-head.fa -o abundances.csv -k 21".split()
    print("** running estimate_query_abundance")
    assert estimate_query_abundance.main(args) == 0

    abunds = open("abundances.csv", "rt").read()


#        assert 'dory-head.fa,1.0,1.05' in abunds


@pytest_utils.in_tempdir
def test_dory_query_by_hashval(location):
    testdata = relative_file("data/dory-k31-hashval-queries.txt")
    shutil.copyfile(testdata, "dory-k31-hashval-queries.txt")

    copy_dory_catlas()

    # index by hashval
    args = "-k 31 dory_k21/bcalm.unitigs.db dory_k21_r1_mh.pickle"
    assert index_cdbg_by_minhash.main(args.split()) == 0

    args = "-k 31 dory_k21 dory_k21_r1 dory_k21_r1_mh.pickle dory-k31-hashval-queries.txt dory_k21_r1_hashval_k31 --contigs-db dory_k21/bcalm.unitigs.db"
    assert query_by_hashval.main(args.split()) == 0
    assert os.path.exists("dory_k21_r1_hashval_k31/hashval_results.csv")


@pytest_utils.in_tempdir
def test_dory_multifasta_query(location):
    copy_dory_head()
    copy_dory_catlas()
    copy_dory_sig()

    # make k-mer search index - FIXTURE
    args = "-k 21 dory_k21_r1 --contigs-db dory_k21/bcalm.unitigs.db".split()
    assert index_cdbg_by_kmer.main(args) == 0

    # index by multifasta
    os.mkdir("dory_k21_r1_multifasta")
    args = "dory_k21 dory_k21_r1 dory_k21_r1_multifasta/multifasta.pickle --query dory-head.fa"
    assert index_cdbg_by_multifasta.main(args.split()) == 0

    args = "-k 21 --scaled 100 dory_k21/bcalm.unitigs.db dory_k21_r1_multifasta/hashval.pickle"
    assert index_cdbg_by_minhash.main(args.split()) == 0

    args = "--hashvals dory_k21_r1_multifasta/hashval.pickle --multi-idx dory_k21_r1_multifasta/multifasta.pickle  --query-sig dory-subset.fq.sig --output dory_k21_r1_multifasta/query-results.csv -k 21 --scaled 100"
    assert query_multifasta_by_sig.main(args.split()) == 0

    args = "--multi-idx dory_k21_r1_multifasta/multifasta.pickle --output-cdbg-record dory_k21_r1_multifasta/multifasta.cdbg_by_record.csv --output-cdbg-annot dory_k21_r1_multifasta/multifasta.cdbg_annot.csv --info-csv dory_k21/contigs.info.csv"
    assert extract_cdbg_by_multifasta.main(args.split()) == 0

    assert os.path.exists("dory_k21_r1_multifasta/multifasta.cdbg_by_record.csv")
    assert os.path.exists("dory_k21_r1_multifasta/multifasta.cdbg_annot.csv")
    assert os.path.exists("dory_k21_r1_multifasta/query-results.csv")


@pytest_utils.in_tempdir
def test_dory_multifasta_annot_x(location):
    copy_dory_head()
    copy_dory_catlas()

    queryfile = relative_file('data/dory-prot-query.faa')

    # index by multifasta
    os.mkdir("dory_k21_r1_multifasta")
    args = f"dory_k21 dory_k21_r1 dory_k21_r1_multifasta/multifasta_x.pickle --query {queryfile} -k 10"
    assert index_cdbg_by_multifasta_x.main(args.split()) == 0


@pytest_utils.in_tempdir
def test_dory_shadow_extract(location):
    # run extract_nodes_by_shadow_ratio
    copy_dory_catlas()

    # make k-mer search index
    args = "dory_k21 dory_k21_r1 shadow_out --contigs-db dory_k21/bcalm.unitigs.db".split()
    print("** running extract_nodes_by_shadow_ratio")
    assert extract_nodes_by_shadow_ratio.main(args) == 0


@pytest_utils.in_tempdir
def test_dory_protein_search(location):
    # run query_by_prot
    copy_dory_catlas()

    queryfile = relative_file('data/dory-prot-query.faa')

    args = f"{queryfile} dory_k21/bcalm.unitigs.db".split()
    query_by_prot.main(args)

    with gzip.open('dory-prot-query.faa.nodes.gz') as fp:
        lines = fp.readlines()
        assert int(lines[0].strip()) == 145
        assert int(lines[1].strip()) == 63
        assert len(lines) == 2


@pytest_utils.in_tempdir
def test_dory_translate_search(location):
    # run query_by_prot with --query-is-dna
    copy_dory_catlas()

    queryfile = relative_file('data/dory-dna-translate-query.fa')

    args = f"{queryfile} dory_k21/bcalm.unitigs.db --query-is-dna".split()
    query_by_prot.main(args)

    with gzip.open('dory-dna-translate-query.fa.nodes.gz') as fp:
        lines = fp.readlines()
        assert int(lines[0].strip()) == 145
        assert int(lines[1].strip()) == 63
        assert len(lines) == 2


@pytest_utils.in_tempdir
def test_extract_neighborhoods_by_cdbg_ids(location):
    # run extract_neighborhoods_by_cdbg_ids
    copy_dory_catlas()

    with gzip.open('node-list.txt.gz', 'wt') as fp:
        print('145', file=fp)
        print('63', file=fp)

    cmdline = 'dory_k21 dory_k21_r1 node-list.txt.gz -o xyz.out'.split()
    extract_neighborhoods_by_cdbg_ids.main(cmdline)

    with gzip.open('xyz.out', 'rt') as fp:
        lines = fp.readlines()
        nodes = set( [int(x.strip()) for x in lines] )

        assert 63 in nodes
        assert 145 in nodes
        assert len(nodes) == 2


@pytest_utils.in_tempdir
def test_dory_cdbg_and_catlas_abund(location):
    # run count_dominator_abundance
    copy_dory_catlas()
    copy_dory_subset()

    args = "dory_k21 dory_k21_r1 dory-subset.fa".split()
    print("** running count_dominator_abundance")
    assert count_dominator_abundance.main(args) == 0

    assert os.path.exists('dory_k21_r1_abund')

    with open('dory_k21_r1_abund/dory-subset.fa.cdbg_abund.csv', newline="") as fp:
        lines = fp.readlines()
        assert len(lines) == 737

        total_abund = 0
        total_node_size = 0
        for line in lines[1:]:
            line = line.strip().split(',')
            total_abund += int(line[1])
            total_node_size += int(line[2])

        assert total_abund == 240920
        assert total_node_size == 212475

    with open('dory_k21_r1_abund/dory-subset.fa.dom_abund.csv', newline="") as fp:
        lines = fp.readlines()
        assert len(lines) == 2836

        for line in lines:
            line = line.strip().split(',')
            if line[1] == '7':             # top dom node
                assert line[2] == '240920' # total k-mer abundances
                assert line[3] == '212475' # total number of k-mers
