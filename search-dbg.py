#! /usr/bin/env python
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

        mxt_dict[node] = mh

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
    p.add_argument('catlas_dir')
    p.add_argument('mh_file')
    p.add_argument('label_list', type=str,
                   help='list of labels that should correspond to MinHash')
    args = p.parse_args()

    basename = os.path.basename(args.catlas_dir)

    original_graph = '%s.gxt' % (basename)
    original_graph = os.path.join(args.catlas_dir, original_graph)
    orig_to_labels = load_orig_to_labels(original_graph)

    graph_mxt = '%s.mxt' % (basename)
    graph_mxt = os.path.join(args.catlas_dir, graph_mxt)
    mxt_dict = load_mxt_dict(graph_mxt)

    query_mh = load_mh_dump(args.mh_file)

    search_labels = set([ int(i) for i in args.label_list.split(',') ])

    ## everything is loaded! collect info & search.

    all_nodes = set(orig_to_labels.keys())
    all_labels = set()
    for vv in orig_to_labels.values():
        all_labels.update(vv)

    print('search labels:', search_labels)
    print('all labels:', all_labels)

    found = set()
    for node_id, subject_mh in mxt_dict.items():
        if query_mh.count_common(subject_mh) > 0:
            found.add(node_id)

    print('found:', len(found))

    pos_nodes = all_nodes.intersection(found)
    neg_nodes = all_nodes - pos_nodes

    print(len(pos_nodes), len(neg_nodes))

    def has_search_label(node_id):
        node_labels = orig_to_labels.get(node_id, set())
        return bool(search_labels.intersection(node_labels))

    # true positives: how many nodes did we find that had the right label?
    tp = 0

    # false positives: how many nodes did we find that didn't have right label?
    fp = 0

    for node_id in pos_nodes:
        if has_search_label(node_id):
            tp += 1
        else:
            fp += 1

    # true negatives: how many nodes did we miss that didn't have right label?
    tn = 0

    # false negatives: how many nodes did we miss that did have right label?
    fn = 0

    for node_id in neg_nodes:
        if not has_search_label(node_id):
            tn += 1
        else:
            fn += 1

    print('')
    print('tp:', tp)
    print('fp:', fp)
    print('fn:', fn)
    print('tn:', tn)
    print('')
    sens = (100.0 * tp / (tp + fn))
    spec = (100.0 * tn / (tn + fp))
    print('sensitivity: %.1f' % (100.0 * tp / (tp + fn)))
    print('specificity: %.1f' % (100.0 * tn / (tn + fp)))

    assert tp + fp + fn + tn == len(all_nodes)

if __name__ == '__main__':
    main()
