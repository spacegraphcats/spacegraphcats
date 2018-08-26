import os.path
import tempfile
import shutil
import screed

from spacegraphcats.catlas import catlas
from spacegraphcats.index import index_contigs_by_kmer
from spacegraphcats.search import extract_nodes_by_query
from spacegraphcats.search import characterize_catlas_regions
from spacegraphcats.search import extract_unassembled_nodes
from spacegraphcats.search import catlas_info
from spacegraphcats.search import extract_contigs
from spacegraphcats.search import estimate_query_abundance
from spacegraphcats.utils import make_bgzf
from spacegraphcats.cdbg import label_cdbg
from spacegraphcats.search import extract_reads
import sourmash_lib


class TempDirectory(object):
    def __init__(self):
        self.tempdir = tempfile.mkdtemp(prefix='sgc_test')

    def __enter__(self):
        return self.tempdir

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            shutil.rmtree(self.tempdir, ignore_errors=True)
        except OSError:
            pass

        if exc_type:
            return False


def relative_filename(filename):
    thisdir = os.path.dirname(__file__)
    pkgdir = os.path.join(thisdir, '../..')
    newpath = os.path.join(pkgdir, filename)
    return os.path.abspath(newpath)


class Args(object):
    pass


def test_dory():
    with TempDirectory() as location:
        from spacegraphcats.cdbg import bcalm_to_gxt

        # make the output directory
        try:
            os.mkdir('dory_k21_r1')
        except FileExistsError:
            pass

        # convert the bcalm file to gxt
        args = ['-k', '-21', '-P',
                relative_filename('dory/bcalm.dory.k21.unitigs.fa'),
                'dory_k21_r1/cdbg.gxt',
                'dory_k21_r1/contigs.fa.gz']

        bcalm_to_gxt.main(args)

        # build catlas

        args = Args()
        args.no_checkpoint = True
        args.level = 0
        args.radius = 1
        args.project = 'dory_k21_r1'

        catlas.main(args)

        # make k-mer search index
        args = '-k 21 dory_k21_r1'.split()
        index_contigs_by_kmer.main(args)

        # do search!!
        args='dory_k21_r1 dory_k21_r1_search_oh0 --query data/dory-head.fa -k 21 --overhead=0.0'.split()

        try:
            extract_nodes_by_query.main(args)
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
            assert last_line == 'data/dory-head.fa,1.0,1.0,1671,2,21,1631,1.0,0.0,0.0'

        # run characterize_catlas_regions
        args = 'dory_k21_r1 dory_k1_r1.vec'.split()
        characterize_catlas_regions.main(args)

        # run extract_unassembled_regions
        args = 'dory_k21_r1 data/dory-head.fa dory.regions -k 21'.split()
        extract_unassembled_nodes.main(args)

        # run catlas info
        catlas_info.main(['dory_k21_r1'])

        # run extract_contigs
        args = ['dory_k21_r1',
                'dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.txt.gz',
                '-o',
                'dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.contigs.fa.gz']
        extract_contigs.main(args)

        assert os.path.exists('dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.contigs.fa.gz')

        # run make_bgzf
        args = ['data/dory-subset.fa', '-o', 'dory/dory.reads.bgz']
        make_bgzf.main(args)

        # run label_cdbg
        args = ['dory_k21_r1',
                'dory/dory.reads.bgz', 'dory_k21_r1/reads.bgz.labels']
        label_cdbg.main(args)

        # run extract_reads
        args = ['dory/dory.reads.bgz',
                'dory_k21_r1/reads.bgz.labels',
                'dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.txt.gz',
                '-o',
                'dory_k21_r1_search_oh0/dory-head.fa.cdbg_ids.reads.fa.gz']
        extract_reads.main(args)

        # calculate query abundances
        args = 'dory_k21_r1 data/dory-head.fa -o abundances.csv -k 21'.split()
        estimate_query_abundance.main(args)

        abunds = open('abundances.csv', 'rt').read()
        assert 'data/dory-head.fa,1.0,1.05' in abunds
