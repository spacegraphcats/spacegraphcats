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

    best_match = -1
    best_match_node = None
    best_mh = None
    for catlas_node, subject_mh in mxt_dict.items():
        match = query_mh.compare(subject_mh)
        if match > best_match:
            best_match = match
            best_match_node = catlas_node
            best_mh = subject_mh

    print('best match: similarity %.3f, catlas node %d' % (best_match,
                                                           best_match_node))
    leaves = catlas.find_level0(best_match_node)

    print('found %d domgraph leaves under catlas node %d' % (len(leaves),
                                                             best_match_node))

    ### finally, count the matches/mismatches between MinHash-found nodes
    ### and expected labels.
    
    orig_nodes = []
    for k in leaves:
        orig_nodes += dom_to_orig[k]

    found_label_counts = {}
    for k in set(orig_nodes):
        for label in orig_to_labels.get(k, []):
            found_label_counts[label] = found_label_counts.get(label, 0) + 1

    all_label_counts = {}
    for k, labels in orig_to_labels.items():
        for label in labels:
            all_label_counts[label] = all_label_counts.get(label, 0) + 1

    label_list = [ int(i) for i in args.label_list.split(',') ]
    all_labels = list(all_label_counts.keys())

    tp = 0
    fp = 0
    fn = 0
    tn = 0

    for k in all_labels:
        if k in label_list:
            tp += found_label_counts.get(k, 0)
            fn += all_label_counts[k] - found_label_counts.get(k, 0)
            
        if k not in label_list:
            fp = found_label_counts.get(k, 0)
            tn += all_label_counts[k] - found_label_counts.get(k, 0)

    print('')
    print('looking for labels:', " ".join([str(i) for i in label_list]))
    print('actually found:', found_label_counts)
    print('all label counts:', all_label_counts)

    print('')

    print('tp:', tp)
    print('fp:', fp)
    print('fn:', fn)
    print('tn:', tn)

    assert tp + fp + fn + tn == sum(all_label_counts.values())

    sys.exit(0)


if __name__ == '__main__':
    main()
