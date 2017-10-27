#! /usr/bin/env python
"""
Some code that Titus would like to play with more but doesn't want to explain
right now.

(Briefly, looks for portions of the graph that have lots of shadow nodes
but very few k-mers; extracts reads from it.)
"""
import argparse
import os
import pickle
import sys
import leveldb
import shutil
import collections
import math

from spacegraphcats.logging import log

import screed
import sourmash_lib
from sourmash_lib import MinHash, signature
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import SigLeaf, search_minhashes
import leveldb
import pickle
from search.frontier_search import (frontier_search, compute_overhead, find_shadow)
import sqlite3
from . import bgzf
#from .search_utils import read_bgzf
from . import search_utils
from .search_utils import get_minhashdb_name


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


def get_minhash_size(node_id, minhash_db):
    mh = load_minhash(node_id, minhash_db)
    if mh:
        return len(mh.get_mins())
    return 0


def main(args=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('--maxsize', type=float, default=1000)
    p.add_argument('readsfile')
    p.add_argument('labeled_reads_sqlite')
    p.add_argument('output')
    p.add_argument('-k', '--ksize', default=31, type=int,
                        help='k-mer size (default: 31)')
    args = p.parse_args(args)

    print('maxsize: {:g}'.format(args.maxsize))

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # load catlas DAG
    top_node_id, dag, dag_up, dag_levels = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    # load mapping between dom nodes and cDBG/graph nodes:
    layer0_to_cdbg = load_layer0_to_cdbg(catlas, domfile)
    print('loaded {} layer 0 catlas nodes'.format(len(layer0_to_cdbg)))

    # load minhash DB
    db_path = get_minhashdb_name(args.catlas_prefix, args.ksize, 0, 0)
    if not db_path:
        print('** ERROR, minhash DB does not exist for {}'.format(args.ksize),
              file=sys.stderr)
        sys.exit(-1)
    minhash_db = leveldb.LevelDB(db_path)

    # calculate the cDBG shadow sizes for each catlas node.
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

    # find highest nodes with shadow size less than given max_size
    def find_terminal_nodes(node_id, path_d, max_size):
        node_list = set()
        for sub_id in dag[node_id]:
            # shadow size
            size = node_shadow_sizes[sub_id]

            if size < max_size:
                node_list.add(sub_id)
            else:
                children = find_terminal_nodes(sub_id, path_d, max_size)
                node_list.update(children)
            
        return node_list

    MAXSIZE = args.maxsize #/ mh.scaled
    counts = collections.defaultdict(int)

    print('finding terminal nodes for {}.'.format(MAXSIZE))
    terminal = find_terminal_nodes(top_node_id, counts, MAXSIZE)
    print('...got {}'.format(len(terminal)))

    # now, go through and collect minhashes for each of them, and
    # estimate total k-mers underneath; filter terminal nodes.
    x = []
    merge_mh = load_minhash(top_node_id, minhash_db).copy_and_clear()
    n_merged = 0
    new_node_set = set()
    for node_id in terminal:
        # get minhash size under this node
        mh = load_minhash(node_id, minhash_db)
        if mh:
            mh_size = len(mh.get_mins()) * mh.scaled
        else:
            mh_size = 1

        # retrieve shadow size, calculate mh_size / shadow_size.
        shadow_size = node_shadow_sizes[node_id]
        ratio = math.log(mh_size, 2) - math.log(shadow_size, 2)

        level = dag_levels[node_id]

        # track basic info
        x.append((ratio, node_id, shadow_size, mh_size, level))

        # merge those with larger shadow than k-mers (somewhat arbitrary :)
        if 2**ratio < 10:
            new_node_set.add(node_id)

            if mh_size > 1:
                n_merged += 1
                merge_mh.merge(mh)

    terminal = new_node_set

    print('terminal node stats for MAXSIZE: {:g}'.format(MAXSIZE))
    print('n merged:', n_merged)
    print('n tnodes:', len(terminal), '; k-mer cover:', len(merge_mh.get_mins()) * merge_mh.scaled)
    mh = load_minhash(top_node_id, minhash_db)
    print('total k-mers:', len(mh.get_mins()) * mh.scaled)

    x.sort(reverse=True)
    for (k, v, a, b, c) in x:
        print('ratio: {:.3f}'.format(2**k), '/ shadow size:', a, '/ kmers:', b, '/ level:', c)

    print("removing empty catlas nodes from the terminal nodes...")
    nonempty_terminal = search_utils.remove_empty_catlas_nodes(terminal,
                                                               minhash_db)
    print("...went from {} to {}".format(len(terminal), len(nonempty_terminal)))
    terminal = nonempty_terminal

    #### extract reads

    cdbg_shadow = set()
    terminal_shadow = find_shadow(terminal, dag)
    for x in terminal_shadow:
        cdbg_shadow.update(layer0_to_cdbg.get(x))

    print('loading graph & labels/foo...')
    
    # sql
    dbfilename = args.labeled_reads_sqlite
    assert os.path.exists(dbfilename), 'sqlite file {} does not exist'.format(dbfilename)
    db = sqlite3.connect(dbfilename)
    cursor = db.cursor()

    total_bp = 0
    watermark_size = 1e7
    watermark = watermark_size
    no_tags = 0
    no_tags_bp = 0
    total_seqs = 0
    output_seqs = 0

    outfp = open(args.output, 'wt')
    reads_grabber = search_utils.GrabBGZF_Random(args.readsfile)

    ## get last offset:
    last_offset = search_utils.sqlite_get_max_offset(cursor)

    print('running query...')
    for n, offset in enumerate(search_utils.sqlite_get_offsets(cursor, cdbg_shadow)):
        if n % 10000 == 0:
            print('...at n {} ({:.1f}% of {})'.format(n, offset / last_offset * 100, args.readsfile), end='\r')

        (record, xx) = reads_grabber.get_sequence_at(offset)
        assert xx == offset, (xx, offset)

        name = record.name
        sequence = record.sequence

        total_bp += len(sequence)
        total_seqs += 1

        outfp.write('>{}\n{}\n'.format(name, sequence))

    print('')
    print('fetched {} reads, {} bp matching terminal.'.format(total_seqs, total_bp))


if __name__ == '__main__':
    main()
