import collections
import csv
import os
import sqlite3

import numpy
from screed.screedRecord import Record
from screed.utils import to_str
from sourmash_lib import MinHash

import bbhash

from spacegraphcats.utils.bgzf.bgzf import BgzfReader


def load_layer1_to_cdbg(cdbg_to_catlas, domfile):
    "Load the mapping between first layer catlas and the original DBG nodes."

    # mapping from cdbg dominators to dominated nodes.
    #domset = {}
    # mapping from catlas node IDs to cdbg nodes
    layer1_to_cdbg = {}


    fp = open(domfile, 'rt')
    for line in fp:
        dom_node, *beneath = line.strip().split(' ')

        dom_node = int(dom_node)
        beneath = map(int, beneath)

        equiv_cdbg_to_catlas = cdbg_to_catlas[dom_node]
        layer1_to_cdbg[equiv_cdbg_to_catlas] = set(beneath)

    fp.close()

    return layer1_to_cdbg


def load_dag(catlas_file):
    "Load the catlas Directed Acyclic Graph."
    dag = {}
    dag_up = collections.defaultdict(set)
    dag_levels = {}
    cdbg_to_catlas = {}

    # track the root of the tree
    max_node = -1
    max_level = -1

    # load everything from the catlas file
    for line in open(catlas_file, 'rt'):
        catlas_node, cdbg_id, level, beneath = line.strip().split(',')

        level = int(level)

        # parse out the children
        catlas_node = int(catlas_node)
        beneath = beneath.strip()
        if beneath:
            beneath = beneath.split(' ')
            beneath = set(map(int, beneath))

            # save node -> children, and level
            dag[catlas_node] = beneath
            for child in beneath:
                dag_up[child].add(catlas_node)
        else:
            dag[catlas_node] = set()

        dag_levels[catlas_node] = level

        # update max_node/max_level
        level = int(level)
        if level > max_level:
            max_level = level
            max_node = catlas_node

        # save cdbg_to_catlas mapping
        if level == 1:
            cdbg_to_catlas[int(cdbg_id)] = catlas_node

    return max_node, dag, dag_up, dag_levels, cdbg_to_catlas


def load_just_dag(catlas_file):
    "Load the catlas Directed Acyclic Graph."
    dag = {}
    cdbg_to_catlas = {}

    # track the root of the tree
    max_node = -1
    max_level = -1

    # load everything from the catlas file
    for line in open(catlas_file, 'rt'):
        catlas_node, cdbg_id, level, beneath = line.strip().split(',')

        level = int(level)

        # parse out the children
        catlas_node = int(catlas_node)
        beneath = beneath.strip()
        if beneath:
            beneath = beneath.split(' ')
            beneath = set(map(int, beneath))

            # save node -> children, and level
            dag[catlas_node] = beneath
        else:
            dag[catlas_node] = set()

        # update max_node/max_level
        level = int(level)
        if level > max_level:
            max_level = level
            max_node = catlas_node

        # save cdbg_to_catlas mapping
        if level == 1:
            cdbg_to_catlas[int(cdbg_id)] = catlas_node

    return max_node, dag, cdbg_to_catlas

# TODO: delete me?
class CatlasDB(object):
    """
    Wrapper class for accessing catlas DAG in sqlite.
    """
    def __init__(self, catlas_prefix):
        self.prefix = catlas_prefix

        dbfilename = 'catlas.csv.sqlite'
        self.dbfilename = os.path.join(self.prefix, dbfilename)

        self.db = sqlite3.connect(self.dbfilename)
        self.cursor = self.db.cursor()
        self.cursor.execute('PRAGMA cache_size=1000000')
        self.cursor.execute('PRAGMA synchronous = OFF')
        self.cursor.execute('PRAGMA journal_mode = MEMORY')
        self.cursor.execute('PRAGMA LOCKING_MODE = EXCLUSIVE')

        self.cursor.execute('BEGIN TRANSACTION')
        self.cursor.execute('SELECT DISTINCT COUNT(node_id) FROM children')
        self.num_nodes = list(self.cursor)[0][0]
        self.cache = collections.defaultdict(set)

    def get_children(self, node_id):
        if node_id in self.cache:
            return self.cache[node_id]
        c = self.cursor
        c.execute('SELECT child_id FROM children WHERE node_id=?', (node_id,))
        xx = set([ x[0] for x in c])
        self.cache[node_id] = xx
        return xx

    def __getitem__(self, node_id):
        x = self.get_children(node_id)
        return x

    def __len__(self):
        return self.num_nodes

    def get_top_node_id(self):
        c = self.cursor
        c.execute('SELECT node_id FROM top_node')
        top_nodes = list(c)
        assert len(top_nodes) == 1
        return top_nodes[0][0]

    def get_shadow(self, node_id_list):
        shadow = set()
        seen_nodes = set()

        def add_to_shadow(node_id):
            if node_id in seen_nodes:
                return
            seen_nodes.add(node_id)

            children_ids = self.get_children(node_id)
            if not len(children_ids):
                shadow.add(node_id)
            else:
                for child in children_ids:
                    add_to_shadow(child)

        for node in node_id_list:
            add_to_shadow(node)

        return shadow

    def get_cdbg_nodes(self, leaf_node_list):
        c = self.cursor()
        cdbg_nodes = set()
        for leaf_node in leaf_node_list:
            c.execute('SELECT cdbg_id FROM first_doms WHERE node_id=?', (leaf_node,))
            for cdbg_id in c:
                cdbg_nodes.update(cdbg_id)

        return cdbg_nodes


def calc_node_shadow_sizes(dag, dag_levels, layer1_to_cdbg):
    x = []
    for (node_id, level) in dag_levels.items():
        x.append((level, node_id))
    x.sort()

    node_shadow_sizes = {}
    for level, node_id in x:
        if level == 0:
            node_shadow_sizes[node_id] = len(layer1_to_cdbg[node_id])
        else:
            sub_size = 0
            for child_id in dag[node_id]:
                sub_size += node_shadow_sizes[child_id]
            node_shadow_sizes[node_id] = sub_size

    return node_shadow_sizes


def sqlite_get_max_offset(cursor):
    cursor.execute('SELECT max(sequences.offset) FROM sequences')
    last_offset, = next(cursor)

    return last_offset


def sqlite_get_offsets(cursor, cdbg_ids):
    seen_offsets = set()
    seen_labels = set()

    cursor.execute('DROP TABLE IF EXISTS label_query')
    cursor.execute('CREATE TEMPORARY TABLE label_query (label_id INTEGER PRIMARY KEY);')

    for label in cdbg_ids:
        cursor.execute('INSERT INTO label_query (label_id) VALUES (?)', (label,))

    cursor.execute('SELECT DISTINCT sequences.offset,sequences.label FROM sequences WHERE label in (SELECT label_id FROM label_query) ORDER BY offset')

    for n, (offset, label) in enumerate(cursor):
        if offset not in seen_offsets:
            yield offset
        seen_offsets.add(offset)
        seen_labels.add(label)

    seen_labels -= cdbg_ids
    assert not seen_labels                # should have gotten ALL the labels


class GrabBGZF_Random(object):
    def __init__(self, filename):
        self.reader = BgzfReader(filename, 'rt')

        ch = self.reader.read(1)

        if ch == '>':
            iter_fn = my_fasta_iter
        elif ch == '@':
            iter_fn = my_fastq_iter
        else:
            raise Exception('unknown start chr {}'.format(ch))

        self.iter_fn = iter_fn

    def get_sequence_at(self, pos):
        self.reader.seek(pos)
        record = next(self.iter_fn(self.reader))
        return record


def iterate_bgzf(reader):
    ch = reader.read(1)

    if ch == '>':
        iter_fn = my_fasta_iter
    elif ch == '@':
        iter_fn = fastq_iter
    else:
        raise Exception('unknown start chr {}'.format(ch))

    reader.seek(0)

    for record, pos in iter_fn(reader):
        yield record, pos


def my_fasta_iter(handle, parse_description=False, line=None):
    """
    Iterator over the given FASTA file handle, returning records. handle
    is a handle to a file opened for reading
    """
    last_start = handle.tell()
    if line is None:
        line = handle.readline()

    while line:
        data = {}

        line = to_str(line.strip())
        if not line.startswith('>'):
            raise IOError("Bad FASTA format: no '>' at beginning of line: {}".format(line))

        if parse_description:  # Try to grab the name and optional description
            try:
                data['name'], data['description'] = line[1:].split(' ', 1)
            except ValueError:  # No optional description
                data['name'] = line[1:]
                data['description'] = ''
        else:
            data['name'] = line[1:]
            data['description'] = ''

        data['name'] = data['name'].strip()
        data['description'] = data['description'].strip()

        # Collect sequence lines into a list
        sequenceList = []
        pos = handle.tell()
        line = to_str(handle.readline())
        while line and not line.startswith('>'):
            sequenceList.append(line.strip())
            pos = handle.tell()
            line = to_str(handle.readline())

        data['sequence'] = ''.join(sequenceList)
        yield Record(**data), last_start
        last_start = pos


def my_fastq_iter(handle, line=None, parse_description=False):
    """
    Iterator over the given FASTQ file handle returning records. handle
    is a handle to a file opened for reading
    """
    if line is None:
        line = handle.readline()
    line = to_str(line.strip())
    while line:
        data = {}

        if line and not line.startswith('@'):
            raise IOError("Bad FASTQ format: no '@' at beginning of line")

        # Try to grab the name and (optional) annotations
        if parse_description:
            try:
                data['name'], data['annotations'] = line[1:].split(' ', 1)
            except ValueError:  # No optional annotations
                data['name'] = line[1:]
                data['annotations'] = ''
                pass
        else:
            data['name'] = line[1:]
            data['annotations'] = ''

        # Extract the sequence lines
        sequence = []
        line = to_str(handle.readline().strip())
        while line and not line.startswith('+') and not line.startswith('#'):
            sequence.append(line)
            line = to_str(handle.readline().strip())

        data['sequence'] = ''.join(sequence)

        # Extract the quality lines
        quality = []
        line = to_str(handle.readline().strip())
        seqlen = len(data['sequence'])
        aclen = 0
        while not line == '' and aclen < seqlen:
            quality.append(line)
            aclen += len(line)
            line = to_str(handle.readline().strip())

        data['quality'] = ''.join(quality)
        if len(data['sequence']) != len(data['quality']):
            raise IOError('sequence and quality strings must be '
                          'of equal length')

        yield Record(**data)


def get_reads_by_cdbg(sqlite_filename, reads_filename, cdbg_ids):
    """
    Given a list of cDBG IDs, retrieve the actual sequences corresponding
    to them by using offsets into a BGZF file.
    """
    # connect to sqlite db
    db = sqlite3.connect(sqlite_filename)
    cursor = db.cursor()

    # open readsfile for random access
    reads_grabber = GrabBGZF_Random(reads_filename)

    ## get last offset in file as measure of progress
    last_offset = sqlite_get_max_offset(cursor)

    # pull out the offsets of all sequences with matches in cdbg_ids.
    for offset in sqlite_get_offsets(cursor, cdbg_ids):
        offset_f = offset / last_offset

        record, xx = reads_grabber.get_sequence_at(offset)
        assert xx == offset

        yield record, offset_f


def get_contigs_by_cdbg(contigs_filename, cdbg_ids):
    """
    Given a list of cDBG IDs, retrieve the actual contig sequences
    corresponding to them by using offsets into a BGZF file.

    This works by iterating over the contig offsets in contigs.fa.gz.info.csv,
    which is created by bcalm_to_gxt; and then seeking into the contigs.fa.gz
    BGZF file to extract the specific sequence.
    """
    info_filename = contigs_filename + '.info.csv'
    reads_grabber = GrabBGZF_Random(contigs_filename)

    with open(info_filename, 'rt') as info_fp:
        r = csv.DictReader(info_fp)

        for row in r:
            contig_id = int(row['contig_id'])
            if contig_id in cdbg_ids:
                offset = int(row['offset'])

                record, xx = reads_grabber.get_sequence_at(offset)
                assert xx == offset
                assert int(record.name) == contig_id, (record.name,contig_id)

                yield record


### MPHF stuff

class MPHF_KmerIndex(object):
    "Support kmer -> cDBG node id queries enabled by search.contigs_by_kmer."
    def __init__(self, mphf, mphf_to_kmer, mphf_to_cdbg_id, cdbg_sizes):
        self.mphf = mphf
        self.mphf_to_kmer = mphf_to_kmer
        self.mphf_to_cdbg_id = mphf_to_cdbg_id
        self.cdbg_sizes = cdbg_sizes

    def get_cdbg_size(self, cdbg_id):
        return self.cdbg_sizes[cdbg_id]

    def get_cdbg_id(self, kmer_hash):
        mphf_hash = self.mphf.lookup(kmer_hash)
        if mphf_hash is None:
            return None
        elif self.mphf_to_kmer[mphf_hash] == kmer_hash:
            cdbg_id = self.mphf_to_cdbg_id[mphf_hash]
            return cdbg_id
        return None

    def build_catlas_node_sizes(self, dag, dag_levels, layer1_to_cdbg):
        x = []
        for (node_id, level) in dag_levels.items():
            x.append((level, node_id))
        x.sort()

        node_kmer_sizes = {}
        for level, node_id in x:
            if level == 1:
                total_kmers = 0
                for cdbg_node in layer1_to_cdbg.get(node_id):
                    total_kmers += self.get_cdbg_size(cdbg_node)

                node_kmer_sizes[node_id] = total_kmers
            else:
                sub_size = 0
                for child_id in dag[node_id]:
                    sub_size += node_kmer_sizes[child_id]
                node_kmer_sizes[node_id] = sub_size

        return node_kmer_sizes

    def get_match_counts(self, query_kmers):
        "Return a dictionary containing cdbg_id -> # of matches in query_kmers"
        match_counts = {}
        n_matched = 0
        for n, hashval in enumerate(query_kmers):
            if n % 1000000 == 0:
                print('matching ...', n, end='\r')

            cdbg_id = self.get_cdbg_id(hashval)
            if cdbg_id is not None:
                match_counts[cdbg_id] = match_counts.get(cdbg_id, 0) + 1
                n_matched += 1

        print('... found {} matches to {} k-mers total.'.format(n_matched,
                                                                n + 1))

        return match_counts

    def build_catlas_match_counts(self, match_counts, dag, dag_levels, layer1_to_cdbg):
        x = []
        for (node_id, level) in dag_levels.items():
            x.append((level, node_id))
        x.sort()

        node_match_counts = {}
        total_kmers_in_query_nodes = 0.
        for level, node_id in x:
            if level == 1:
                query_kmers = 0
                all_kmers_in_node = 0
                for cdbg_node in layer1_to_cdbg.get(node_id):
                    query_kmers += match_counts.get(cdbg_node, 0)
                    all_kmers_in_node += self.get_cdbg_size(cdbg_node)

                if query_kmers:
                    node_match_counts[node_id] = query_kmers
                    total_kmers_in_query_nodes += all_kmers_in_node

                    assert all_kmers_in_node >= query_kmers

            else:
                sub_size = 0
                for child_id in dag[node_id]:
                    sub_size += node_match_counts.get(child_id, 0)

                if sub_size:
                    node_match_counts[node_id] = sub_size

        return node_match_counts


def load_kmer_index(catlas_prefix):
    "Load kmer index created by search.contigs_by_kmer."
    mphf_filename = os.path.join(catlas_prefix, 'contigs.fa.gz.mphf')
    array_filename = os.path.join(catlas_prefix, 'contigs.fa.gz.indices')
    mphf = bbhash.load_mphf(mphf_filename)
    with open(array_filename, 'rb') as fp:
        np_dict = numpy.load(fp)

        mphf_to_kmer = np_dict['mphf_to_kmer']
        mphf_to_cdbg = np_dict['kmer_to_cdbg']
        cdbg_sizes = np_dict['sizes']

    return MPHF_KmerIndex(mphf, mphf_to_kmer, mphf_to_cdbg, cdbg_sizes)


def load_cdbg_size_info(catlas_prefix, min_abund=0.0):
    filename = os.path.join(catlas_prefix, 'contigs.fa.gz.info.csv')
    with open(filename, 'rt') as fp:
        cdbg_kmer_sizes = {}
        cdbg_weighted_kmer_sizes = {}
        r = csv.DictReader(fp)
        for row in r:
            contig_id = int(row['contig_id'])
            n_kmers = int(row['n_kmers'])
            mean_abund = float(row['mean_abund'])
            if not min_abund or mean_abund >= min_abund:
                cdbg_kmer_sizes[contig_id] = n_kmers
                cdbg_weighted_kmer_sizes[contig_id] = mean_abund*n_kmers

    return cdbg_kmer_sizes, cdbg_weighted_kmer_sizes


def decorate_catlas_with_shadow_sizes(layer1_to_cdbg, dag, dag_levels):
    x = []
    for (node_id, level) in dag_levels.items():
        x.append((level, node_id))
    x.sort()

    node_shadow_sizes = {}
    for level, node_id in x:
        if level == 1:
            node_shadow_sizes[node_id] = len(layer1_to_cdbg[node_id])
        else:
            sub_size = 0
            for child_id in dag[node_id]:
                sub_size += node_shadow_sizes[child_id]
            node_shadow_sizes[node_id] = sub_size

    return node_shadow_sizes


def decorate_catlas_with_kmer_sizes(layer1_to_cdbg, dag, dag_levels, cdbg_kmer_sizes, cdbg_weighted_kmer_sizes):
    x = []
    for (node_id, level) in dag_levels.items():
        x.append((level, node_id))
    x.sort()

    node_kmer_sizes = {}
    node_weighted_kmer_sizes = {}
    for level, node_id in x:
        if level == 1:                    # aggregate across cDBG nodes
            total_kmers = 0
            total_weighted_kmers = 0
            for cdbg_node in layer1_to_cdbg.get(node_id):
                total_kmers += cdbg_kmer_sizes.get(cdbg_node, 0)
                total_weighted_kmers += cdbg_weighted_kmer_sizes.get(cdbg_node, 0)

            node_kmer_sizes[node_id] = total_kmers
            node_weighted_kmer_sizes[node_id] = total_weighted_kmers
        else:                             # aggregate across children
            sub_size = 0
            sub_weighted_size = 0
            for child_id in dag[node_id]:
                sub_size += node_kmer_sizes[child_id]
                sub_weighted_size += node_weighted_kmer_sizes[child_id]
            node_kmer_sizes[node_id] = sub_size
            node_weighted_kmer_sizes[node_id] = sub_weighted_size

    return node_kmer_sizes, node_weighted_kmer_sizes


def output_response_curve(outname, match_counts, kmer_idx, layer1_to_cdbg):
    curve = []

    # track total containment
    total = 0

    # walk over all layer1 nodes
    for node_id, cdbg_nodes in sorted(layer1_to_cdbg.items()):
        n_matches = 0
        n_kmers = 0

        # aggregate counts across cDBG nodes under this layer1 node.
        for cdbg_node in cdbg_nodes:
            n_matches += match_counts.get(cdbg_node, 0)
            n_kmers += kmer_idx.get_cdbg_size(cdbg_node)

        # do we keep this layer1 node, i.e. does it have positive containment?
        if n_matches:
            n_cont = n_matches
            n_oh = (n_kmers - n_matches)

            total += n_cont
            curve.append((n_cont, n_oh, node_id))

    # sort by absolute containment
    curve.sort(reverse=True)

    # track total containment etc
    sofar = 0
    total_oh = 0
    total_cont = 0

    # @CTB: remove redundant sum_cont fields
    # @CTB: switch to CSV output
    # @CTB: ask Mike what he wants here :)

    with open(outname, 'wt') as fp:
        fp.write('sum_cont relative_cont relative_overhead sum_cont2 sum_oh catlas_id\n')

        # only output ~200 points
        sampling_rate = max(int(len(curve) / 200), 1)

        # do the output thing
        for pos, (n_cont, n_oh, node_id) in enumerate(curve):
            sofar += n_cont
            total_oh += n_oh
            total_cont += n_cont

            # output first and last points, as well as at sampling rate.
            if pos % sampling_rate == 0 or pos == 0 or pos+1 == len(curve):
                fp.write('{} {} {} {} {} {}\n'.format(sofar, total_cont / total, total_oh / total, n_cont, n_oh, node_id))
