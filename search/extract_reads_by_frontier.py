#! /usr/bin/env python
import argparse
import os
import sys
import leveldb
from copy import copy
import collections

import sourmash_lib
from sourmash_lib import MinHash, signature
from sourmash_lib.sourmash_args import load_query_signature
from typing import Dict, List, Set, Union, Tuple
from pickle import load
import screed

from .memoize import memoize
from .search_catlas_with_minhash import load_dag, load_minhash
from spacegraphcats.logging import log
from search.frontier_search import (frontier_search, compute_overhead, find_shadow)
from .search_utils import (load_dag, load_layer0_to_cdbg)
import khmer.utils


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query_sig', help='query minhash')
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('overhead', help='\% of overhead', type=float)
    p.add_argument('labeled_reads')
    p.add_argument('labeled_reads_labels')
    p.add_argument('output')
    p.add_argument('--no-empty', action='store_true')
    p.add_argument('--purgatory', action='store_true')
    p.add_argument('--fullstats', action='store_true')
    p.add_argument('-k', '--ksize', default=None, type=int,
                        help='k-mer size (default: 31)')

    args = p.parse_args()

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # load catlas DAG
    top_node_id, dag, dag_up, dag_levels = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    db_path = os.path.join(args.catlas_prefix, 'minhashes.db')
    minhash_db = leveldb.LevelDB(db_path)

    # load mapping between dom nodes and cDBG/graph nodes:
    layer0_to_cdbg = load_layer0_to_cdbg(catlas, domfile)
    print('loaded {} layer 0 catlas nodes'.format(len(layer0_to_cdbg)))
    x = set()
    for v in layer0_to_cdbg.values():
        x.update(v)
    print('...corresponding to {} cDBG nodes.'.format(len(x)))

    # load query MinHash
    query_sig = load_query_signature(args.query_sig, select_ksize=args.ksize,
                                     select_moltype='DNA')
    print('loaded query sig {}'.format(query_sig.name()))

    frontier, num_leaves, num_empty, frontier_mh = frontier_search(query_sig, top_node_id, dag, minhash_db, args.overhead, not args.no_empty, args.purgatory)

    top_mh = load_minhash(top_node_id, minhash_db)
    query_mh = query_sig.minhash.downsample_max_hash(top_mh)
    top_mh = top_mh.downsample_max_hash(query_sig.minhash)
    print("Root containment: {}".format(query_mh.contained_by(top_mh)))
    print("Root similarity: {}".format(query_mh.similarity(top_mh)))

    print("Containment of frontier: {}".format(query_mh.contained_by(frontier_mh)))
    print("Similarity of frontier: {}".format(query_mh.similarity(frontier_mh)))
    print("Size of frontier: {} of {} ({:.3}%)".format(len(frontier), len(dag), 100 * len(frontier) / len(dag)))
    print("Overhead of frontier: {}".format(compute_overhead(frontier_mh, query_mh)))
    print("Number of leaves in the frontier: {}".format(num_leaves))
    print("Number of empty catlas nodes in the frontier: {}".format(num_empty))
    print("")

    print("removing...")
    nonempty_frontier = []
    merge_mh = query_mh.copy_and_clear()
    max_scaled = max(merge_mh.scaled, top_mh.scaled)

    for node in frontier:
        mh = load_minhash(node, minhash_db)
        if mh and len(mh.get_mins()) > 0:
            nonempty_frontier.append(node)
            mh = mh.downsample_scaled(max_scaled)
            merge_mh.merge(mh)
    print("...went from {} to {}".format(len(frontier), len(nonempty_frontier)))
    print('recalculated frontier mh similarity: {}'.format(merge_mh.similarity(query_mh)))
    frontier = nonempty_frontier

    shadow = find_shadow(frontier, dag)

    print("Size of the frontier shadow: {}".format(len(shadow)))
    if len(shadow) == len(layer0_to_cdbg):
        print('\n*** WARNING: shadow is the entire graph! ***\n')

    query_size = len(query_sig.minhash.get_mins())
    query_bp = query_size * query_sig.minhash.scaled
    print("Size of query minhash: {} (est {:2.1e} bp)".\
              format(query_size, query_bp))
    minhash_size = len(frontier_mh.get_mins())
    minhash_bp = minhash_size * frontier_mh.scaled
    print("Size of frontier minhash: {} (est {:2.1e} bp); ratio {:.2f}".\
              format(minhash_size, minhash_bp, minhash_bp / query_bp))

    log(args.catlas_prefix, sys.argv)

    #### extract reads

    cdbg_shadow = set()
    for x in shadow:
        cdbg_shadow.update(layer0_to_cdbg.get(x))

    print('loading graph & labels/foo...')
    # actually don't need whole graph, just tags @CTB
    ng = khmer.Nodegraph(args.ksize, 1, 1)
    ng.load_tagset(args.labeled_reads_labels + '.tagset')
    tags_to_labels = load(open(args.labeled_reads_labels + '.labelsp', 'rb'))

    total_bp = 0
    watermark_size = 1e7
    watermark = watermark_size
    no_tags = 0
    no_tags_bp = 0
    total_seqs = 0
    output_seqs = 0

    outfp = open(args.output, 'wt')

    n = 0
    n_no_tags = 0
    n_no_labels = 0
    empty_set = set()
    for record in screed.open(args.labeled_reads):
        n += 1
        if total_bp >= watermark:
            print('... {:5.2e} bp thru reads; {} read, {} written, {} no tags'.format(int(watermark), total_seqs, output_seqs, n_no_tags),
                  file=sys.stderr, end='\r')
            watermark += watermark_size

        total_bp += len(record.sequence)
        total_seqs += 1

        tags = ng.get_tags_for_sequence(record.sequence)

        labels = set()
        if tags:
            for t in tags:
                labels.update(tags_to_labels.get(t, empty_set))
        else:
            n_no_tags += 1

        if not labels:
            n_no_labels += 1
            continue

        if labels.intersection(cdbg_shadow):
            khmer.utils.write_record(record, outfp)
            output_seqs += 1

    print('')
    print(no_tags, total_seqs)
    print('{:5.2e} {:5.2e}'.format(no_tags_bp, total_bp))
    

    sys.exit(0)


if __name__ == '__main__':
    # import cProfile
    # cProfile.run('main()', 'search_stats')

    main()
