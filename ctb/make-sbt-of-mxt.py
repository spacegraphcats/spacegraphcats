#! /usr/bin/env python
import os
import sys
import argparse
from sourmash_lib import MinHash
from collections import defaultdict
from spacegraphcats.catlas import CAtlas
import time
import sourmash_lib
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import search_minhashes, SigLeaf
from sourmash_lib import signature


MINHASH_K=31
SCALED=100.0

def import_graph_mxt(mxt):
    "Load all of the minhashes from an MXT file into a dict."
    d = defaultdict(list)
    with open(mxt, 'rt') as fp:
        for line in fp:
            line = line.strip().split(' ')
            g_id = int(line[0])
            if len(line) > 1:
                mins = line[1]

                mins = list(map(int, mins.split(' ')))
                d[g_id] = mins

    return d


def load_layer1_to_cdbg(catlas_file):
    "Load the mapping between first layer catlas and the original DBG nodes."
    layer1_to_cdbg = defaultdict(set)

    catlas_to_cdbg = {}
    for line in open(catlas_file, 'rt'):
        catlas_node, cdbg_node, level, beneath = line.strip().split(',')
        if int(level) != 0:
            continue

        catlas_node = int(catlas_node)
        cdbg_node = int(cdbg_node)
        catlas_to_cdbg[catlas_node] = cdbg_node

    for line in open(catlas_file, 'rt'):
        catlas_node, cdbg_node, level, beneath = line.strip().split(',')
        if int(level) != 1:
            continue

        catlas_node = int(catlas_node)
        beneath = beneath.split(' ')
        beneath = list(map(int, beneath))

        for catnode in beneath:
            dbgnode = catlas_to_cdbg[catnode]
            layer1_to_cdbg[catlas_node].add(dbgnode)

    return layer1_to_cdbg


def merge_nodes(child_dict, child_node_list):
    """Merge child nodes into a single minhash."""
    minlist = []

    for graph_node in child_node_list:
        if graph_node in child_dict:
            mins = child_dict[graph_node]
            minlist.extend(mins)

    # merge into a single minhash!
    merged_mh = MinHash(0, MINHASH_K,
                        max_hash=round(sourmash_lib.MAX_HASH / SCALED))
    for h in minlist:
        merged_mh.add_hash(h)

    assert len(merged_mh.get_mins()) == len(set(minlist))

    # add into merged minhashes table.
    return merged_mh


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('-x', '--bf-size', type=float, default=1e4)
    args = p.parse_args()
    
    basename = os.path.basename(args.catlas_prefix)
    graphmxt = '%s.mxt' % (basename,)
    graphmxt = os.path.join(args.catlas_prefix, graphmxt)

    catlas = os.path.join(args.catlas_prefix, args.catlas_prefix + '.catlas')
    
    # load MXT into dict
    graph_minhashes = import_graph_mxt(graphmxt)
    print('imported {} graph minhashes'.format(len(graph_minhashes)))

    # load mapping between dom nodes and cDBG/graph nodes:
    layer1_to_cdbg = load_layer1_to_cdbg(catlas)
    print('loaded {} layer 1 catlas nodes'.format(len(layer1_to_cdbg)))
    x = set()
    for v in layer1_to_cdbg.values():
        x.update(v)
    print('...corresponding to {} cDBG nodes.'.format(len(x)))

    #print(layer1_to_cdbg[1100])
    #assert(1199 in layer1_to_cdbg[1100])
    #assert(1201 in layer1_to_cdbg[1100])
    #assert(1514 in layer1_to_cdbg[1100])

    # create minhashes for catlas leaf nodes.
    leaf_minhashes = {}
    for n, (catlas_node, cdbg_nodes) in enumerate(layer1_to_cdbg.items()):
        mh = merge_nodes(graph_minhashes, cdbg_nodes)
        leaf_minhashes[catlas_node] = mh
    print('created {} leaf node MinHashes via merging'.format(n + 1))

    factory = GraphFactory(1, args.bf_size, 4)
    tree = SBT(factory)

    print('building tree...')
    for node_id, mh in leaf_minhashes.items():
        ss = signature.SourmashSignature('', mh, name='{}'.format(node_id))

        leaf = SigLeaf(ss.md5sum(), ss)
        tree.add_node(leaf)

    print('...done with {} minhashes. saving!'.format(len(leaf_minhashes)))
    sbt_name = os.path.basename(args.catlas_prefix)
    tree.save(sbt_name)
    print('saved sbt "{}"'.format(sbt_name))

    xxx_mh = MinHash(0, MINHASH_K,
                        max_hash=round(sourmash_lib.MAX_HASH / SCALED))
    n = 0
    for minlist in graph_minhashes.values():
        for hash in minlist:
            xxx_mh.add_hash(hash)
        n += 1
    print('zzz', n, len(graph_minhashes))

    xx = signature.SourmashSignature('', xxx_mh, name='all')
    with open('all.sig', 'wt') as fp:
        signature.save_signatures([xx], fp)

    sys.exit(0)

if __name__ == '__main__':
    main()
