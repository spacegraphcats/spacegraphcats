#! /usr/bin/env python
"""
Make the minhashes for catlas nodes, based on the cDBG nodes they dominate.
"""
import argparse
import os
import pickle
import sys
import shutil
from collections import defaultdict
import gc

from spacegraphcats.logging import log
from search import search_utils

import screed
from sourmash_lib import MinHash


class MinHashFactory(object):
    def __init__(self, **params):
        self.params = params

    def __call__(self):
        return MinHash(**self.params)


def make_leaf_minhashes(contigfile, cdbg_to_layer1, factory):
    "Make the minhashes for each leaf node from the contigs in contigfile"

    d = defaultdict(factory)
    total_bp = 0
    watermark = 1e7

    fp = screed.open(contigfile)
    for record in fp:
        if total_bp >= watermark:
            print('... {:5.2e} bp thru contigs'.format(int(watermark)),
                  file=sys.stderr)
            watermark += 1e7

        # the FASTA header/record name is the cDBG node ID output by
        # build_contracted_dbg.
        cdbg_id = int(record.name)

        # make a minhash!
        mh = factory()
        mh.add_sequence(record.sequence)
        if mh.get_mins():
            leaf_nodes = cdbg_to_layer1.get(cdbg_id, set())
            for leaf_node_id in leaf_nodes:
                d[leaf_node_id].merge(mh)

        total_bp += len(record.sequence)

    fp.close()

    return d


def build_catlas_minhashes(catlas_file, catlas_minhashes, factory, save_db):
    "Build MinHashes for all the internal nodes of the catlas DAG."

    total_mh = 0
    empty_mh = 0

    # create a list of all the nodes, sorted by level (increasing)
    x = []
    fp = open(catlas_file, 'rt')
    for line in fp:
        catlas_node, cdbg_node, level, beneath = line.strip().split(',')
        if int(level) == 1:
            continue

        beneath = beneath.split(' ')
        beneath = list(map(int, beneath))

        x.append((int(level), int(catlas_node), beneath))
    fp.close()

    # walk through, building the merged minhashes (which we can do in a
    # single pass on the sorted list).
    x.sort()
    # In general, the level of v's children is one less than v's level.  This
    # is not true, however if v is the root, since we add vertices to be
    # children of the root as soon as they become isolated in the domgraph.
    # Our cleanup of no-longer needed minhashes can't delete the children of
    # the root until all other nodes have been processed.
    root_children = set(x[-1][2])

    levels = defaultdict(set)
    current_level = 1
    for (level, catlas_node, beneath) in x:

        # remove no-longer needed catlas minhashes to save on memory
        if level > current_level:
            if level > 2:
                print('flushing catlas minhashes at level {}'.format(level-2))
                for node in levels[level - 2]:
                    # don't delete yet if it's a child of the root
                    if node not in root_children:
                        del catlas_minhashes[node]

            current_level = level

        # track catlas nodes thus far for later deletion
        levels[level].add(catlas_node)

        # merge!
        merged_mh = merge_nodes(catlas_minhashes, beneath, factory)

        if not merged_mh.get_mins():
            merged_mh = None

        catlas_minhashes[catlas_node] = merged_mh
        total_mh += 1

        # write!
        if merged_mh:
            if save_db:
                save_db.save_mh(catlas_node, merged_mh)
        else:
            empty_mh += 1

    return total_mh, empty_mh


def merge_nodes(child_dict, child_node_list, factory):
    """Merge child nodes into a single minhash."""
    # merge into a single minhash!
    merged_mh = factory()

    for graph_node in child_node_list:
        if graph_node in child_dict:
            mh = child_dict[graph_node]
            if mh:
                mh = mh.downsample_scaled(merged_mh.scaled)
                merged_mh.merge(mh)

    # add into merged minhashes table.
    return merged_mh


def main(args=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('--track-abundance', action='store_true')
    p.add_argument('--seed', default=42,
                   type=int)
    p.add_argument('--num', default=10, type=int)

    args = p.parse_args(args)

    ksize = args.ksize
    seed = args.seed
    track_abundance = args.track_abundance

    # build factories to produce new MinHash objects.
    factory = MinHashFactory(n=args.num, ksize=ksize,
                             track_abundance=args.track_abundance,
                             seed=args.seed)
    leaf_factory = factory

    # put together the basic catlas info --
    basename = os.path.basename(args.catlas_prefix)
    contigfile = os.path.join(args.catlas_prefix, "contigs.fa.gz")

    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # load mapping between dom nodes and cDBG/graph nodes:
    max_node, dag, cdbg_to_catlas = search_utils.load_just_dag(catlas)
    del dag
    gc.collect()

    layer1_to_cdbg = search_utils.load_layer1_to_cdbg(cdbg_to_catlas, domfile)
    print('loaded {} layer 1 catlas nodes'.format(len(layer1_to_cdbg)))
    x = set()
    for v in layer1_to_cdbg.values():
        x.update(v)
    print('...corresponding to {} cDBG nodes.'.format(len(x)))
    del x

    # build reverse mapping
    cdbg_to_layer1 = defaultdict(set)
    for catlas_node, cdbg_nodes in layer1_to_cdbg.items():
        for cdbg_id in cdbg_nodes:
            cdbg_to_layer1[cdbg_id].add(catlas_node)

    # create the minhash db, first removing it if it already exists
    path = os.path.join(args.catlas_prefix, 'minhash_FOO.db')

    print('saving minhashes in {}'.format(path))
    if os.path.exists(path):
        print('(removing previously existing file)')
        os.unlink(path)
    save_db = search_utils.MinhashSqlDB(path)
    save_db.create()

    # create minhashes for catlas leaf nodes.
    print('ksize={} num={}'.format(ksize, args.num))
    catlas_minhashes = make_leaf_minhashes(contigfile, cdbg_to_layer1,
                                           leaf_factory)

    total_mh = len(layer1_to_cdbg)
    empty_mh = total_mh - len(catlas_minhashes)
    print('... built {} leaf node MinHashes.'.format(total_mh),
          file=sys.stderr)

    for catlas_node, mh in catlas_minhashes.items():
        save_db.save_mh(catlas_node, mh)

    # build minhashes for entire catlas
    t, e = build_catlas_minhashes(catlas, catlas_minhashes, factory,
                                  save_db)
    total_mh += t
    empty_mh += e

    print('committing DB...')
    save_db.commit()

    print('saved {} minhashes (and {} empty)'.format(total_mh - empty_mh,
                                                     empty_mh))

    # log that this command was run
    log(args.catlas_prefix, sys.argv)


if __name__ == '__main__':
    main()
