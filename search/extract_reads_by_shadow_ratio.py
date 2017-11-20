#! /usr/bin/env python
"""
Some code that Titus would like to play with more but doesn't want to explain
right now.

(Briefly, looks for portions of the graph that have lots of shadow nodes
but very few k-mers; extracts reads from it.)
"""
import argparse
import os
import sys
import leveldb
import collections
import math

from . import search_utils
from .search_utils import get_minhashdb_name, get_reads_by_cdbg, load_minhash
from .frontier_search import find_shadow


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
    top_node_id, dag, dag_up, dag_levels, cdbg_to_catlas = search_utils.load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    # load mapping between dom nodes and cDBG/graph nodes:
    layer1_to_cdbg = search_utils.load_layer1_to_cdbg(cdbg_to_catlas, domfile)
    print('loaded {} layer 1 catlas nodes'.format(len(layer1_to_cdbg)))

    # load minhash DB
    db_path = get_minhashdb_name(args.catlas_prefix, args.ksize, 0, 0, 43)
    if not db_path:
        print('** ERROR, minhash DB does not exist for {}'.format(args.ksize),
              file=sys.stderr)
        sys.exit(-1)
    minhash_db = search_utils.MinhashSqlDB(db_path)

    # calculate the cDBG shadow sizes for each catlas node.
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
        print(shadow_size)
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

    # build cDBG shadow ID list.
    cdbg_shadow = set()
    terminal_shadow = find_shadow(terminal, dag)
    for x in terminal_shadow:
        cdbg_shadow.update(layer1_to_cdbg.get(x))

    #### extract reads
    print('loading graph & labels/foo...')
    
    dbfilename = args.labeled_reads_sqlite
    assert os.path.exists(dbfilename), 'sqlite file {} does not exist'.format(dbfilename)

    total_bp = 0
    total_seqs = 0
    output_seqs = 0

    outfp = open(args.output, 'wt')

    print('running query...')
    reads_iter = get_reads_by_cdbg(dbfilename, args.readsfile, cdbg_shadow)
    for n, (record, offset_f) in enumerate(reads_iter):
        if n % 10000 == 0:
            print('...at n {} ({:.1f}% of {})'.format(n, offset_f * 100,
                                                      args.readsfile),
                  end='\r')

        total_bp += len(record.sequence)
        total_seqs += 1

        outfp.write('>{}\n{}\n'.format(record.name, record.sequence))

    print('')
    print('fetched {} reads, {} bp matching terminal.'.format(total_seqs, total_bp))


if __name__ == '__main__':
    main()
