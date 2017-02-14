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

def import_graph_mxt(mxt):
    "Load all of the minhashes from an MXT file into a dict."
    d = {}
    with open(mxt, 'rt') as fp:
        for line in fp:
            line = line.split(',')
            g_id = int(line[0])
            mins = line[1]

            mins = list(map(int, mins.split(' ')))
            d[g_id] = mins

    return d


def load_dom_to_orig(assignment_vxt_file):
    "Load the mapping between dom nodes and the original DBG nodes."
    dom_to_orig = defaultdict(list)
    for line in open(assignment_vxt_file, 'rt'):
        orig_node, dom_list = line.strip().split(',')
        orig_node = int(orig_node)
        dom_list = list(map(int, dom_list.split(' ')))

        for dom_node in dom_list:
            dom_to_orig[dom_node].append(orig_node)

    return dom_to_orig


def merge_nodes(child_dict, child_node_list):
    """Merge child nodes into a single minhash."""
    minlist = []

    for graph_node in child_node_list:
        mins = child_dict[graph_node]
        minlist.append(mins)

    # merge into a single minhash!
    merged_mh = MinHash(len(minlist), MINHASH_K)
    for mins in minlist:
        for h in mins:
            merged_mh.add_hash(h)

    assert len(merged_mh) == len(minlist)

    # add into merged minhashes table.
    return merged_mh


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('catlas_r', help='catlas radius', type=int)
    p.add_argument('-x', '--bf-size', type=float, default=1e4)
    args = p.parse_args()
    
    radius = args.catlas_r
    basename = os.path.basename(args.catlas_prefix)
    graphmxt = '%s.mxt' % (basename,)
    graphmxt = os.path.join(args.catlas_prefix, graphmxt)
    
    catgxt = '%s.catlas.%d.gxt' % (basename, radius)
    catmxt = '%s.catlas.%d.mxt' % (basename, radius)
    catgxt = os.path.join(args.catlas_prefix, catgxt)
    catmxt = os.path.join(args.catlas_prefix, catmxt)

    catmxt_db = catmxt + '.db'

    # load MXT into dict
    graph_minhashes = import_graph_mxt(graphmxt)
    print('imported {} graph minhashes'.format(len(graph_minhashes)))

    # load mapping between dom nodes and cDBG/graph nodes:
    assignment_vxt = '%s.assignment.%d.vxt' % (basename, radius)
    assignment_vxt = os.path.join(args.catlas_prefix, assignment_vxt)
    dom_to_orig = load_dom_to_orig(assignment_vxt)

    # create minhashes for catlas leaf nodes.
    dom_minhashes = {}
    for n, (dom_node, graph_nodes) in enumerate(dom_to_orig.items()):
        mh = merge_nodes(graph_minhashes, graph_nodes)
        dom_minhashes[dom_node] = mh
    print('created {} leaf node MinHashes via merging'.format(n + 1))

    factory = GraphFactory(1, args.bf_size, 4)
    tree = SBT(factory)

    print('building tree...')
    for node_id, mh in dom_minhashes.items():
        e = sourmash_lib.Estimators(ksize=MINHASH_K, n=len(mh))
        e.mh = mh
        ss = signature.SourmashSignature('', e, name='{}'.format(node_id))

        leaf = SigLeaf(ss.md5sum(), ss)
        tree.add_node(leaf)

    print('...done with {} minhashes. saving!'.format(len(dom_minhashes)))
    sbt_name = os.path.basename(args.catlas_prefix) + '.' + str(radius)
    tree.save(sbt_name)
    print('saved sbt "{}"'.format(sbt_name))

    sys.exit(0)

if __name__ == '__main__':
    main()
