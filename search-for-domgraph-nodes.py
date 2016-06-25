#! /usr/bin/env python3
"""
"""
from __future__ import print_function
import argparse
from khmer import MinHash
from spacegraphcats import graph_parser
from spacegraphcats.catlas_reader import CAtlasReader
import os
import sys


KSIZE=31


def load_orig_to_labels(original_graph_filename):
    "Load the labels for the original DBG node IDs"
    orig_to_labels = {}
    def parse_source_graph_labels(node_id, size, names, vals):
        assert names[0] == 'labels'
        labels = vals[0]
        if labels:
            orig_to_labels[node_id] = list(map(int, labels.strip().split(' ')))
        
    def nop(*x):
        pass

    with open(original_graph_filename) as fp:
        graph_parser.parse(fp, parse_source_graph_labels, nop)

    return orig_to_labels


def load_dom_to_orig(assignment_vxt_file):
    "Load the mapping between level0/dom nodes and the original DBG nodes."
    dom_to_orig = {}
    for line in open(assignment_vxt_file, 'rt'):
        orig_node, dom_list = line.strip().split(',')
        orig_node = int(orig_node)
        dom_list = list(map(int, dom_list.split(' ')))

        for dom_node in dom_list:
            x = dom_to_orig.get(dom_node, [])
            x.append(orig_node)
            dom_to_orig[dom_node] = x

    return dom_to_orig


def load_mxt_dict(mxt_filename):
    "Load the MXT file containing minhashes for each DAG node."
    mxt_dict = {}
    for line in open(mxt_filename):
        node, hashes = line.strip().split(',')
        node = int(node)
        hashes = [ int(h) for h in hashes.split(' ') ]
        mh = MinHash(len(hashes), KSIZE)
        for h in hashes:
            mh.add_hash(h)
            
        mxt_dict[int(node)] = mh

    return mxt_dict


def load_mh_dump(mh_filename):
    "Load a MinHash from a file created with 'sourmash dump'."
    with open(mh_filename) as fp:
        hashes = fp.read().strip().split(' ')
    query_mh = MinHash(len(hashes), KSIZE)

    for h in hashes:
        query_mh.add_hash(int(h))
    return query_mh


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('catlas_r', type=int, help='catlas radius to load')
    p.add_argument('mh_file', help='file containing dumped MinHash signature')
    p.add_argument('label_list', type=str,
                   help='list of labels that should correspond to MinHash')
    args = p.parse_args()

    ### first, parse the catlas gxt

    catlas = CAtlasReader(args.catlas_prefix, args.catlas_r)

    ### get the labels from the original graph

    # here, 'orig_to_labels' is a dictionary mapping De Bruijn graph node IDs
    # to a list of labels for each node.
    #
    #   orig_to_labels[dbg_node_id] => list of [label_ids]

    orig_to_labels = load_orig_to_labels(catlas.original_graph)

    ### backtrack the leaf nodes to the domgraph

    # here, 'dom_to_orig' is a dictionary mapping domination nodes to
    # a list of original vertices in the De Bruijn graph.
    #
    #   dom_to_orig[dom_node_id] => list of [orig_node_ids]

    dom_to_orig = load_dom_to_orig(catlas.assignment_vxt)

    ### load mxt

    # 'mxt_dict' is a dictionary mapping catlas node IDs to MinHash
    # objects.

    print('reading mxt file', catlas.catlas_mxt)
    mxt_dict = load_mxt_dict(catlas.catlas_mxt)

    ### load mh

    print('reading mh file', args.mh_file)

    query_mh = load_mh_dump(args.mh_file)

    ### next, find the relevant catlas nodes using the MinHash.

    print('searching catlas minhashes w/%s' % args.mh_file)

    match_nodes = catlas.find_matching_nodes_best_match(query_mh, mxt_dict)

    leaves = set()
    for match_node in match_nodes:
        leaves.update(catlas.find_level0(match_node))

    print('found %d domgraph leaves under catlas nodes %s' % \
              (len(leaves), match_nodes))

    ### finally, count the matches/mismatches between MinHash-found nodes
    ### and expected labels.

    search_labels = set([ int(i) for i in args.label_list.split(',') ])
    all_labels = set()
    for vv in orig_to_labels.values():
        all_labels.update(vv)

    print("labels we're looking for:", search_labels)
    print("all labels:", all_labels)

    # all_nodes is set of all labeled node_ids on original graph:
    all_nodes = set(dom_to_orig.keys())

    # pos_nodes is set of MH-matching node_ids from original cDBG.
    pos_nodes = set()
    for k in leaves:
        pos_nodes.update(dom_to_orig[k])

    # intersect pos_nodes with all_nodes so we only have labeled
    pos_nodes.intersection_update(all_nodes)

    # neg_nodes is set of non-MH-matching node IDs
    neg_nodes = all_nodes - pos_nodes

    print('')
    print('pos nodes:', len(pos_nodes))
    print('neg nodes:', len(neg_nodes))
    print('all nodes:', len(all_nodes))

    def has_search_label(node_id):
        node_labels = orig_to_labels.get(node_id, set())
        return bool(search_labels.intersection(node_labels))

    # true positives: how many nodes did we find that had the right label?
    tp = 0

    # false positives: how many nodes did we find that didn't have right label?
    fp = 0
    
    for orig_node_id in pos_nodes:
        if has_search_label(orig_node_id):
            tp += 1
        else:
            fp += 1

    # true negatives: how many nodes did we miss that didn't have right label?
    tn = 0
    
    # false negatives: how many nodes did we miss that did have right label?
    fn = 0
    
    for orig_node_id in neg_nodes:
        if not has_search_label(orig_node_id):
            tn += 1
        else:
            fn += 1

    print('')
    print('tp:', tp)
    print('fp:', fp)
    print('fn:', fn)
    print('tn:', tn)

    assert tp + fp + fn + tn == len(all_nodes)

    sys.exit(0)


if __name__ == '__main__':
    main()
