import pickle
import collections
import os
import json

import leveldb
import sqlite3

import screed
import sourmash_lib
from sourmash_lib import MinHash
from screed.screedRecord import Record
from screed.utils import to_str

from .bgzf.bgzf import BgzfReader


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

    # track the root of the tree
    max_node = -1
    max_level = -1

    # load everything from the catlas file
    for line in open(catlas_file, 'rt'):
        catlas_node, cdbg_node, level, beneath = line.strip().split(',')

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

    return max_node, dag, dag_up, dag_levels


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


def load_minhash(node_id: int, minhash_db: leveldb.LevelDB) -> MinHash:
    "Load an individual node's MinHash from the leveldb."
    try:
        value = minhash_db.Get(node_id.to_bytes(8, byteorder='big'))
    except KeyError:
        return None

    return pickle.loads(value)


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


def remove_empty_catlas_nodes(nodes, minhash_db):
    import time
    start = time.time()
    nonempty_frontier = set()
    for node in nodes:
        mh = load_minhash(node, minhash_db)
        if mh and len(mh.get_mins()) > 0:
            nonempty_frontier.add(node)

    print('removed {} empty catlas nodes ({:.1f} s)'.format(len(nodes) - len(nonempty_frontier), time.time() - start))

    return nonempty_frontier


def boost_frontier(frontier, frontier_mh, dag, dag_up, minhash_db, top_node_id):
    """
    Find an internal frontier that covers the given frontier with no add'l
    overhead.
    """
    boosted_frontier = set()

    for node in frontier:
        while 1:
            node_up = dag_up.get(node)
            if not node_up:               # top!
                assert node == top_node_id, node
                break

            assert len(node_up) == 1          # with minor construction
            node_up = node_up.pop()

            node_up_mh = load_minhash(node_up, minhash_db)
            node_up_mh = node_up_mh.downsample_scaled(frontier_mh.scaled)
            if node_up_mh.contained_by(frontier_mh) < 1.0:
                break

            node = node_up

        # try one more...
        node_up = dag_up[node]
        if node_up:
            node_up = node_up.pop()
            if node_up != top_node_id:
                node_up_mh = load_minhash(node_up, minhash_db)
                node_up_mh = node_up_mh.downsample_scaled(frontier_mh.scaled)
                print('FOO!!', node_up_mh.contained_by(frontier_mh))
                if node_up_mh.contained_by(frontier_mh) > 0.8:
                    node = node_up
            else:
                print('next node up is top!')
        else:
            print('WTF?', node, node_up, top_node_id)

        boosted_frontier.add(node)

    return boosted_frontier


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


def get_minhashdb_name(catlas_prefix, ksize, scaled, track_abundance, seed,
                       must_exist=True):
    """
    Construct / return the name of the minhash db given the parameters.
    """
    track_abundance = bool(track_abundance)

    # first, check if it's in minhashes_info.json
    if must_exist:
        infopath = os.path.join(catlas_prefix, 'minhashes_info.json')
        if not os.path.exists(infopath):
            return None

        info = []
        with open(infopath, 'rt') as fp:
            info = json.loads(fp.read())

        matches = []
        for d in info:
            if d['ksize'] == ksize and d['seed'] == seed and \
              d['track_abundance'] == track_abundance:
                matches.append(d['scaled'])

        if not matches:
            return None

        matches.sort(reverse=True)
        found = False
        for db_scaled in matches:
            if not scaled or db_scaled <= scaled:
                 found = True
                 scaled = db_scaled
                 break

        if not found:
            return None

    # ok, now create name.
    is_abund = 0
    if track_abundance:
        is_abund = 1

    name = 'minhashes.db.k{}.s{}.abund{}.seed{}'
    name = name.format(ksize, scaled, is_abund, seed)
    path = os.path.join(catlas_prefix, name)

    if must_exist and not os.path.exists(path):
        return None

    return path


def update_minhash_info(catlas_prefix, ksize, scaled, track_abundance, seed):
    """
    Update minhashes_info with new db info.
    """
    infopath = os.path.join(catlas_prefix, 'minhashes_info.json')
    info = []
    if os.path.exists(infopath):
        with open(infopath, 'rt') as fp:
            info = json.loads(fp.read())

    this_info = dict(ksize=ksize, scaled=scaled,
                     track_abundance=track_abundance, seed=seed)
    if this_info not in info:
        info.append(this_info)

        with open(infopath, 'wt') as fp:
            fp.write(json.dumps(info))


def build_queries_for_seeds(seeds, ksize, scaled, query_seq_file):
    seed_mh_list = []

    for seed in seeds:
        mh = MinHash(0, ksize, scaled=scaled, seed=seed)
        seed_mh_list.append(mh)

    name = None
    for record in screed.open(query_seq_file):
        if not name:
            name = record.name
        for seed_mh in seed_mh_list:
            seed_mh.add_sequence(record.sequence, True)

    seed_queries = [sourmash_lib.SourmashSignature(seed_mh, name=name) for \
                        seed_mh in seed_mh_list]
    return seed_queries


def parse_seeds_arg(seeds_str):
    seeds = []
    seeds_str = seeds_str.split(',')
    for seed in seeds_str:
        if '-' in seed:
            (start, end) = seed.split('-')
            for s in range(int(start), int(end) + 1):
                seeds.append(s)
        else:
            seeds.append(int(seed))

    return seeds
