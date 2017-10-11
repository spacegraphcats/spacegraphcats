import pickle
import collections

import leveldb
import sqlite3

from sourmash_lib import MinHash
from screed.screedRecord import Record
from screed.utils import to_str

from .bgzf.bgzf import BgzfReader


def load_layer0_to_cdbg(catlas_file, domfile):
    "Load the mapping between first layer catlas and the original DBG nodes."

    # mapping from cdbg dominators to dominated nodes.
    domset = {}

    fp = open(domfile, 'rt')
    for line in fp:
        dom_node, *beneath = line.strip().split(' ')

        dom_node = int(dom_node)
        beneath = map(int, beneath)

        domset[dom_node] = set(beneath)

    fp.close()

    layer0_to_cdbg = {}

    # mapping from catlas node IDs to cdbg nodes
    fp = open(catlas_file, 'rt')
    for line in fp:
        catlas_node, cdbg_node, level, beneath = line.strip().split(',')
        if int(level) != 0:
            continue

        catlas_node = int(catlas_node)
        cdbg_node = int(cdbg_node)
        layer0_to_cdbg[catlas_node] = domset[cdbg_node]

    fp.close()

    return layer0_to_cdbg


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


def load_minhash(node_id: int, minhash_db: leveldb.LevelDB) -> MinHash:
    "Load an individual node's MinHash from the leveldb."
    try:
        value = minhash_db.Get(node_id.to_bytes(8, byteorder='big'))
    except KeyError:
        return None

    return pickle.loads(value)


def calc_node_shadow_sizes(dag, dag_levels, layer0_to_cdbg):
    x = []
    for (node_id, level) in dag_levels.items():
        x.append((level, node_id))
    x.sort()

    node_shadow_sizes = {}
    for level, node_id in x:
        if level == 0:
            node_shadow_sizes[node_id] = len(layer0_to_cdbg[node_id])
        else:
            sub_size = 0
            for child_id in dag[node_id]:
                sub_size += node_shadow_sizes[child_id]
            node_shadow_sizes[node_id] = sub_size

    return node_shadow_sizes


def remove_empty_catlas_nodes(nodes, minhash_db):
    nonempty_frontier = set()
    for node in nodes:
        mh = load_minhash(node, minhash_db)
        if mh and len(mh.get_mins()) > 0:
            nonempty_frontier.add(node)

    return nonempty_frontier


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
