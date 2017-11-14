#! /usr/bin/env python
"""
Make the minhashes for catlas nodes, based on the cDBG nodes they dominate.
"""
import argparse
import os
import pickle
import sys
import leveldb
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


class LevelDBWriter(object):
    def __init__(self, path):
        self.db = leveldb.LevelDB(path)
        self.batch = None

    def start(self):
        self.batch = leveldb.WriteBatch()

    def end(self):
        assert self.batch
        self.db.Write(self.batch, sync=True)
        self.batch = None

    def put_minhash(self, node_id, mh):
        b = node_id.to_bytes(8, byteorder='big')
        p = pickle.dumps(mh)
        self.batch.Put(b, p)


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
    levels = defaultdict(set)
    current_level = 1
    for (level, catlas_node, beneath) in x:

        # remove no-longer needed catlas minhashes to save on memory
        if level > current_level:
            if level > 2:
                print('flushing catlas minhashes at level {}'.format(level-2))
                for node in levels[level - 2]:
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
                save_db.put_minhash(catlas_node, merged_mh)
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
    p.add_argument('--leaves-only', action='store_true')
    p.add_argument('--scaled', default=1000, type=float)
    p.add_argument('--leaf-scaled', default=1000, type=float)
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('--track-abundance', action='store_true')
    p.add_argument('--seed', default=42,
                   type=int)

    args = p.parse_args(args)

    ksize = args.ksize
    scaled = int(args.scaled)
    seed = args.seed
    leaf_scaled = int(args.leaf_scaled)
    track_abundance = args.track_abundance

    # build factories to produce new MinHash objects.
    factory = MinHashFactory(n=0, ksize=ksize,
                             scaled=scaled,
                             track_abundance=args.track_abundance,
                             seed=args.seed)
    leaf_factory = MinHashFactory(n=0, ksize=ksize,
                             scaled=leaf_scaled,
                             track_abundance=args.track_abundance,
                             seed=args.seed)

    # put together the basic catlas info --
    basename = os.path.basename(args.catlas_prefix)
    contigfile = os.path.join(args.catlas_prefix, "contigs.fa.gz")

    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # load mapping between dom nodes and cDBG/graph nodes:
    layer1_to_cdbg = search_utils.load_layer1_to_cdbg(catlas, domfile)
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
    path = search_utils.get_minhashdb_name(args.catlas_prefix, ksize, scaled,
                                           track_abundance, seed,
                                           must_exist=False)
    if os.path.exists(path):
        shutil.rmtree(path)

    print('saving minhashes in {}'.format(path))
    save_db = LevelDBWriter(path)
    save_db.start()                       # batch mode writing

    # create minhashes for catlas leaf nodes.
    print('ksize={} scaled={}'.format(ksize, scaled))
    catlas_minhashes = make_leaf_minhashes(contigfile, cdbg_to_layer1,
                                           leaf_factory)

    total_mh = len(layer1_to_cdbg)
    empty_mh = total_mh - len(catlas_minhashes)
    print('... built {} leaf node MinHashes.'.format(total_mh),
          file=sys.stderr)

    for catlas_node, mh in catlas_minhashes.items():
        if save_db:
            save_db.put_minhash(catlas_node, mh)

    # build minhashes for entire catlas, or just the leaves (dom nodes)?
    if not args.leaves_only:
        t, e = build_catlas_minhashes(catlas, catlas_minhashes, factory,
                                      save_db)
        total_mh += t
        empty_mh += e

    save_db.end()
    print('saved {} minhashes (and {} empty)'.format(total_mh - empty_mh,
                                                     empty_mh))

    # write out some metadata
    search_utils.update_minhash_info(args.catlas_prefix, ksize, scaled,
                                     track_abundance, seed)

    # log that this command was run
    log(args.catlas_prefix, sys.argv)


if __name__ == '__main__':
    main()
