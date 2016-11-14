#! /usr/bin/env python3
import os
import sys
import argparse
from sourmash_lib import MinHash
from collections import defaultdict
from spacegraphcats.catlas import CAtlas


MINHASH_SIZE=1000
MINHASH_K=0

def import_graph_mxt(mxt):
    "Load all of the minhashes from an MXT file into a dict."
    n = 0

    d = {}
    with open(mxt, 'rt') as fp:
        for line in fp:
            line = line.split(',')
            g_id = int(line[0])
            mins = line[1]

            d[g_id] = list(map(int, mins.split(' ')))

    return n, d


def export_catlas_mxt(catlas_minhashes, fp):
    "Dump all of the minhashes from 'catlas_minhashes' to an MXT file."
    n = 0

    for (node_id, mins) in catlas_minhashes.items():
        fp.write('{},{}\n'.format(node_id, " ".join(map(str, mins))))
        n += 1
        
    return n


def merge_nodes(child_dict, child_node_list):
    """Merge child_node_list from 'from_tablename' into parent_node in
    'to_tablename'."""
    minlist = []

    for graph_node in child_node_list:
        mins = child_dict[graph_node]
        minlist.append(mins)

    # merge!
    # (note that minhashes are implicitly merged when you simply add all
    # the hashes to one MH).
    merged_mh = MinHash(MINHASH_SIZE, MINHASH_K)
    for mins in minlist:
        for h in mins:
            merged_mh.add_hash(h)

    # add into merged minhashes table.
    return merged_mh.get_mins()


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


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('catlas_r', help='catlas radius', type=int)
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

    # load mapping between dom nodes and cDBG/graph nodes:
    assignment_vxt = '%s.assignment.%d.vxt' % (basename, radius)
    assignment_vxt = os.path.join(args.catlas_prefix, assignment_vxt)
    dom_to_orig = load_dom_to_orig(assignment_vxt)
    
    # load MXT into dict
    n, graph_minhashes = import_graph_mxt(graphmxt)
    print('imported {} graph minhashes'.format(n))

    # create minhashes for catlas leaf nodes.
    dom_minhashes = {}
    for n, (dom_node, graph_nodes) in enumerate(dom_to_orig.items()):
        mins = merge_nodes(graph_minhashes, graph_nodes)
        dom_minhashes[dom_node] = mins
    print('created {} leaf node MinHashes via merging'.format(n + 1))

    # this gives us the leaf level minhashes that we need for the rest.
    # now, eliminate the assignment dict & go for the catlas structure.
    del dom_to_orig
    catlas = CAtlas.read(catgxt, None, args.catlas_r)

    # for level 0, merge the shadows (domgraph nodes)
    print('merging level 0')
    n = 0
    m = 0
    select = lambda node: node.level == 0
    catlas_minhashes = {}
    for node in catlas.nodes(select):
        assert node.id not in catlas_minhashes
        domgraph_nodes = node.shadow()
        mins = merge_nodes(dom_minhashes, domgraph_nodes)
        catlas_minhashes[node.id] = mins
        n += 1
        m += len(domgraph_nodes)

    print('level 0: merged {} children into {} nodes'.format(m, n))

    # for each level above 0, merge the children
    for level in range(1, catlas.level + 1):
        print('merging at level:', level)
        n = 0
        m = 0
        select = lambda node: node.level == level
        for node in catlas.nodes(select):
            assert node.id not in catlas_minhashes
            mins = merge_nodes(catlas_minhashes,
                               [ x.id for x in node.children ])
            catlas_minhashes[node.id] = mins
            n += 1
            m += len(node.children)
        print('level {}: merged {} children into {} nodes'.format(level, m, n))

    with open(catmxt, 'wt') as fp:
        n = export_catlas_mxt(catlas_minhashes, fp)
        print('exported {} catlas minhashes to {}'.format(n, catmxt))

    sys.exit(0)

if __name__ == '__main__':
    main()
