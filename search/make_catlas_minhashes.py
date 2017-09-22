#! /usr/bin/env python
import argparse
import os
import pickle
import sys
import leveldb
import shutil
import json
from collections import defaultdict

from spacegraphcats.logging import log

import screed
import sourmash_lib
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


def make_contig_minhashes(contigfile, factory):
    "Make the minhashes for each contig in the contigfile."

    d = {}
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
        if not mh.get_mins():
            mh = None

        d[cdbg_id] = mh

        total_bp += len(record.sequence)

    fp.close()

    return d


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


def build_catlas_minhashes(catlas_file, catlas_minhashes, factory, save_db):
    "Build MinHashes for all the internal nodes of the catlas DAG."

    total_mh = 0
    empty_mh = 0

    # create a list of all the nodes, sorted by level (increasing)
    x = []
    fp = open(catlas_file, 'rt')
    for line in fp:
        catlas_node, cdbg_node, level, beneath = line.strip().split(',')
        if int(level) == 0:
            continue

        beneath = beneath.split(' ')
        beneath = list(map(int, beneath))

        x.append((int(level), int(catlas_node), beneath))
    fp.close()

    # walk through, building the merged minhashes (which we can do in a
    # single pass on the sorted list).
    x.sort()
    levels = defaultdict(set)
    current_level = 0
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
                merged_mh.merge(mh)

    # add into merged minhashes table.
    return merged_mh


def main(args=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('--leaves-only', action='store_true')
    p.add_argument('--scaled', default=100.0, type=float)
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('--track-abundance', action='store_true')

    args = p.parse_args(args)

    ksize = args.ksize
    scaled = args.scaled
    track_abundance = args.track_abundance

    # build a factory to produce new MinHash objects.
    max_hash = sourmash_lib.MAX_HASH / float(scaled)
    max_hash = int(round(max_hash, 0))

    factory = MinHashFactory(n=0, ksize=ksize,
                             max_hash=max_hash,
                             track_abundance=args.track_abundance)

    # put together the basic catlas info --
    basename = os.path.basename(args.catlas_prefix)
    contigfile = os.path.join(args.catlas_prefix, "contigs.fa.gz")

    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # make minhashes from node contigs
    print('ksize={} scaled={:.0f}'.format(ksize, scaled))
    print('making contig minhashes...')
    graph_minhashes = make_contig_minhashes(contigfile, factory)
    print('...made {} contig minhashes'.format(len(graph_minhashes)))

    # load mapping between dom nodes and cDBG/graph nodes:
    layer0_to_cdbg = load_layer0_to_cdbg(catlas, domfile)
    print('loaded {} layer 0 catlas nodes'.format(len(layer0_to_cdbg)))
    x = set()
    for v in layer0_to_cdbg.values():
        x.update(v)
    print('...corresponding to {} cDBG nodes.'.format(len(x)))

    # create the minhash db
    path = os.path.join(args.catlas_prefix, 'minhashes.db')

    if os.path.exists(path):
        shutil.rmtree(path)

    print('saving minhashes in {}'.format(path))
    save_db = LevelDBWriter(path)
    save_db.start()                       # batch mode writing

    # create minhashes for catlas leaf nodes.
    catlas_minhashes = {}
    total_mh = 0
    empty_mh = 0
    for n, (catlas_node, cdbg_nodes) in enumerate(layer0_to_cdbg.items()):
        if n and n % 10000 == 0:
            print('... built {} leaf node MinHashes...'.format(n),
                  file=sys.stderr)
        mh = merge_nodes(graph_minhashes, cdbg_nodes, factory)
        catlas_minhashes[catlas_node] = mh

        total_mh += 1
        if mh:
            if save_db:
                save_db.put_minhash(catlas_node, mh)
        else:
            empty_mh += 1                 # track empty

    print('created {} leaf node MinHashes via merging'.format(n + 1))
    print('')

    # build minhashes for entire catlas, or just the leaves (dom nodes)?
    if not args.leaves_only:
        t, e = build_catlas_minhashes(catlas, catlas_minhashes, factory,
                                      save_db)
        total_mh += t
        empty_mh += e

    save_db.end()
    print('saved {} minhashes ({} empty)'.format(total_mh - empty_mh,
                                                 empty_mh))

    # write out some metadata
    infopath = os.path.join(args.catlas_prefix, 'minhashes_info.json')
    with open(infopath, 'w') as fp:
        info = [ dict(ksize=ksize, scaled=scaled,
                      track_abundance=track_abundance) ]
        fp.write(json.dumps(info))

    # log that this command was run
    log(args.catlas_prefix, sys.argv)

if __name__ == '__main__':
    main()
