import os.path
import shutil
import screed
import sourmash

import pytest
import spacegraphcats.utils.pytest_utils as pytest_utils
from spacegraphcats.utils.pytest_utils import (pkg_file, relative_file)

from spacegraphcats.catlas import catlas
from spacegraphcats.cdbg import index_contigs_by_kmer
from spacegraphcats.search import query_by_sequence
from spacegraphcats.search import characterize_catlas_regions
from spacegraphcats.search import extract_unassembled_nodes
from spacegraphcats.search import evaluate_overhead
from spacegraphcats.search import catlas_info
from spacegraphcats.search import extract_contigs
from spacegraphcats.search import estimate_query_abundance
from spacegraphcats.search import extract_nodes_by_shadow_ratio
from spacegraphcats.utils import make_bgzf
from spacegraphcats.cdbg import label_cdbg
from spacegraphcats.search import extract_reads
from spacegraphcats.cdbg import index_cdbg_by_minhash
from spacegraphcats.search import query_by_hashval
from spacegraphcats.search import index_cdbg_by_multifasta
from spacegraphcats.cdbg import index_cdbg_by_minhash
from spacegraphcats.search import query_multifasta_by_sig
from spacegraphcats.search import extract_cdbg_by_multifasta

def copy_dory_catlas():
    testdata = pkg_file('search/test-data/catlas.dory_k21_r1')
    shutil.copytree(testdata, './dory_k21_r1')

def copy_dory_catlas_search():
    testdata = pkg_file('search/test-data/catlas.dory_k21_r1_search_oh0')
    shutil.copytree(testdata, './dory_k21_r1_search_oh0')

def copy_dory_head():
    testdata = relative_file('data/dory-head.fa')
    shutil.copyfile(testdata, 'dory-head.fa')

def copy_dory_subset():
    testdata = relative_file('data/dory-subset.fa')
    shutil.copyfile(testdata, 'dory-subset.fa')

    testdata = relative_file('data/dory-subset.fq')
    shutil.copyfile(testdata, 'dory-subset.fq')

def copy_dory_sig():
    testdata = relative_file('data/dory-subset.fq.sig')
    shutil.copyfile(testdata, 'dory-subset.fq.sig')

### actual tests

@pytest_utils.in_tempdir
def test_dory_query_workflow(location):
    from spacegraphcats.cdbg import bcalm_to_gxt
    copy_dory_head()
    copy_dory_subset()

    # make the output directory
    try:
        os.mkdir('dory_k21_r1')
    except FileExistsError:
        pass

    # convert the bcalm file to gxt
    args = ['-k', '21', '-P',
            relative_file('dory/bcalm.dory.k21.unitigs.fa'),
            'dory_k21_r1/cdbg.gxt',
            'dory_k21_r1/contigs.fa.gz']

    bcalm_to_gxt.main(args)

    # build catlas
    args = pytest_utils.Args()
    args.no_checkpoint = True
    args.level = 0
    args.radius = 1
    args.project = 'dory_k21_r1'
    print('** running catlas')
    catlas.main(args)

    # make k-mer search index
    args = '-k 21 dory_k21_r1'.split()
    print('** running index_contigs_by_kmer')
    index_contigs_by_kmer.main(args)

    # do search!!
    args='dory_k21_r1 dory_k21_r1_search_oh0 --query dory-head.fa -k 21'.split()
    try:
        query_by_sequence.main(args)
    except SystemExit as e:
        assert e.code == 0, str(e)

    # check output!
    output_path = 'dory_k21_r1_search_oh0/'
    assert os.path.exists(output_path + 'command.txt')
    assert os.path.exists(output_path + 'dory-head.fa.frontier.txt.gz')
    assert os.path.exists(output_path + 'dory-head.fa.cdbg_ids.txt.gz')
    assert os.path.exists(output_path + 'dory-head.fa.response.txt')
    assert os.path.exists(output_path + 'dory-head.fa.contigs.sig')
    assert os.path.exists(output_path + 'results.csv')

    with open(output_path + 'results.csv') as fp:
        lines = fp.readlines()
        assert len(lines) == 2

        last_line = lines[-1].strip()
        assert last_line == 'dory-head.fa,1.0,1.0,1671,2,21,1631,1.0,0.0,0.0,dory_k21_r1'

    sigfile = output_path + 'dory-head.fa.contigs.sig'
    with open(sigfile, 'rt') as fp:
        sig = sourmash.load_one_signature(fp)
        name = sig.name()
        assert 'from dory_k21_r1' in name
        assert name.startswith('nbhd:TRINITY_DN290219_c0_g1_i1')


@pytest_utils.in_tempdir
def test_dory_search_nomatch(location):
    # test situations where zero k-mers match - should not fail.
    copy_dory_catlas()

    testdata = relative_file('data/random-query-nomatch.fa')
    shutil.copyfile(testdata, 'random-query.fa')

    # make k-mer search index
    args = '-k 21 dory_k21_r1'.split()
    print('** running index_contigs_by_kmer')
    index_contigs_by_kmer.main(args)

    # do search!!
    args='dory_k21_r1 dory_k21_r1_search_oh0 --query random-query.fa -k 21'.split()
    try:
        query_by_sequence.main(args)
    except SystemExit as e:
        assert e.code == 0, str(e)


@pytest_utils.in_tempdir
def test_dory_characterize_catlas_regions(location):
    copy_dory_catlas()
    copy_dory_head()
    copy_dory_subset()

    # run characterize_catlas_regions
    args = 'dory_k21_r1 dory_k1_r1.vec'.split()
    characterize_catlas_regions.main(args)


@pytest_utils.in_tempdir
def test_dory_extract_unassembled_nodes(location):
    copy_dory_catlas()
    copy_dory_head()

    # make k-mer search index - FIXTURE
    args = '-k 21 dory_k21_r1'.split()
    index_contigs_by_kmer.main(args)

    # run extract_unassembled_regions
    args = 'dory_k21_r1 dory-head.fa dory.regions -k 21'.split()
    extract_unassembled_nodes.main(args)


@pytest_utils.in_tempdir
def test_dory_catlas_info(location):
    copy_dory_catlas()

    # run catlas info
    print('running catlas_info')
    catlas_info.main(['dory_k21_r1'])


@pytest_utils.in_tempdir
def test_dory_extract_contigs(location):
    copy_dory_catlas()
    copy_dory_catlas_search()
    copy_dory_head()

    # run extract_contigs
    print('running extract_info')
    args = ['dory_k21_r1',
            'dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.txt.gz',
            '-o',
            'dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.contigs.fa.gz']
    extract_contigs.main(args)

    assert os.path.exists('dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.contigs.fa.gz')


@pytest_utils.in_tempdir
def test_dory_make_bgzf(location):
    copy_dory_subset()

    # run make_bgzf
    print('** running make_bgzf')
    args = ['dory-subset.fa', '-o', relative_file('dory/dory.reads.bgz')]
    make_bgzf.main(args)


@pytest_utils.in_tempdir
def test_dory_label_cdbg(location):
    copy_dory_catlas()
    copy_dory_subset()

    # run make_bgzf - FIXTURE
    print('** running make_bgzf')
    args = ['dory-subset.fa', '-o', relative_file('dory/dory.reads.bgz')]
    make_bgzf.main(args)

    # run label_cdbg
    print('** running label_cdbg')
    args = ['dory_k21_r1',
            relative_file('dory/dory.reads.bgz'),
            'dory_k21_r1/reads.bgz.labels']
    label_cdbg.main(args)


@pytest_utils.in_tempdir
def test_dory_extract_reads(location):
    copy_dory_catlas()
    copy_dory_catlas_search()
    copy_dory_subset()

    # run make_bgzf - FIXTURE
    print('** running make_bgzf')
    args = ['dory-subset.fa', '-o', relative_file('dory/dory.reads.bgz')]
    make_bgzf.main(args)

    # run label_cdbg - FIXTURE
    print('** running label_cdbg')
    args = ['dory_k21_r1',
            relative_file('dory/dory.reads.bgz'),
            'dory_k21_r1/reads.bgz.labels']
    label_cdbg.main(args)

    # run extract_reads
    print('** running extract_reads')
    args = [relative_file('dory/dory.reads.bgz'),
            'dory_k21_r1/reads.bgz.labels',
            'dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.txt.gz',
            '-o',
            'dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.reads.fa.gz']
    extract_reads.main(args)

    reads_filename = 'dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.reads.fa.gz'
    reads = [ record for record in screed.open(reads_filename) ]
    assert len(reads) == 2


@pytest_utils.in_tempdir
def test_dory_extract_reads_fq(location):
    # check that FASTQ is passed through properly.

    copy_dory_catlas()
    copy_dory_catlas_search()
    copy_dory_subset()

    # run make_bgzf - FIXTURE
    print('** running make_bgzf')
    args = ['dory-subset.fq', '-o', relative_file('dory/dory.reads.bgz')]
    make_bgzf.main(args)

    # run label_cdbg - FIXTURE
    print('** running label_cdbg')
    args = ['dory_k21_r1',
            relative_file('dory/dory.reads.bgz'),
            'dory_k21_r1/reads.bgz.labels']
    label_cdbg.main(args)

    # run extract_reads
    print('** running extract_reads')
    args = [relative_file('dory/dory.reads.bgz'),
            'dory_k21_r1/reads.bgz.labels',
            'dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.txt.gz',
            '-o',
            'dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.reads.fa.gz']
    extract_reads.main(args)

    reads_filename = 'dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.reads.fa.gz'
    reads = [ record for record in screed.open(reads_filename) ]
    assert len(reads) == 2
    assert len(reads[0].quality)          # FASTQ preserved!


@pytest_utils.in_tempdir
def test_dory_evaluate_overhead(location):
    copy_dory_catlas()
    copy_dory_head()
    copy_dory_catlas_search()

    # run evaluate_overhead
    args = ['dory_k21_r1',
            'dory-head.fa',
            'dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.txt.gz',
            '-o', 'xyz']
    print('** running evaluate_overhead')
    evaluate_overhead.main(args)


@pytest_utils.in_tempdir
def test_dory_estimate_query_abundance(location):
    copy_dory_catlas()
    copy_dory_head()

    # make k-mer search index - FIXTURE
    args = '-k 21 dory_k21_r1'.split()
    index_contigs_by_kmer.main(args)

    # calculate query abundances
    args = 'dory_k21_r1 dory-head.fa -o abundances.csv -k 21'.split()
    print('** running estimate_query_abundance')
    estimate_query_abundance.main(args)

    abunds = open('abundances.csv', 'rt').read()
#        assert 'dory-head.fa,1.0,1.05' in abunds


@pytest_utils.in_tempdir
def test_dory_query_by_hashval(location):
    testdata = relative_file('data/dory-k31-hashval-queries.txt')
    shutil.copyfile(testdata, 'dory-k31-hashval-queries.txt')

    copy_dory_catlas()

    # index by hashval
    args = '-k 31 dory_k21_r1/contigs.fa.gz dory_k21_r1_mh.pickle'
    index_cdbg_by_minhash.main(args.split())

    args = '-k 31 dory_k21_r1 dory_k21_r1_mh.pickle dory-k31-hashval-queries.txt dory_k21_r1_hashval_k31'
    query_by_hashval.main(args.split())
    assert os.path.exists('dory_k21_r1_hashval_k31/hashval_results.csv')


@pytest_utils.in_tempdir
def test_dory_multifasta_query(location):
    copy_dory_head()
    copy_dory_catlas()
    copy_dory_sig()

    # make k-mer search index - FIXTURE
    args = '-k 21 dory_k21_r1'.split()
    index_contigs_by_kmer.main(args)

    # index by multifasta
    os.mkdir('dory_k21_r1_multifasta')
    args = 'dory_k21_r1 dory_k21_r1_multifasta/multifasta.pickle --query dory-head.fa -k 21'
    index_cdbg_by_multifasta.main(args.split())

    args = '-k 21 --scaled 100 dory_k21_r1/contigs.fa.gz dory_k21_r1_multifasta/hashval.pickle'
    index_cdbg_by_minhash.main(args.split())

    args = '--hashvals dory_k21_r1_multifasta/hashval.pickle --multi-idx dory_k21_r1_multifasta/multifasta.pickle  --query-sig dory-subset.fq.sig --output dory_k21_r1_multifasta/query-results.csv -k 21 --scaled 100'
    query_multifasta_by_sig.main(args.split())

    args = '--multi-idx dory_k21_r1_multifasta/multifasta.pickle --output dory_k21_r1_multifasta/multifasta.cdbg_by_record.csv --info-csv dory_k21_r1/contigs.fa.gz.info.csv'
    extract_cdbg_by_multifasta.main(args.split())

    assert os.path.exists('dory_k21_r1_multifasta/multifasta.cdbg_by_record.csv')
    assert os.path.exists('dory_k21_r1_multifasta/query-results.csv')


@pytest_utils.in_tempdir
def test_dory_shadow_extract(location):
    # run extract_nodes_by_shadow_ratio
    copy_dory_catlas()

    # make k-mer search index
    args = 'dory_k21_r1 shadow_out'.split()
    print('** running extract_nodes_by_shadow_ratio')
    extract_nodes_by_shadow_ratio.main(args)
