#! /usr/bin/env python
from __future__ import print_function
import argparse
from sourmash_lib import MinHash
from spacegraphcats import graph_parser
from spacegraphcats.catlas_reader import CAtlasReader
import os
import sys

KSIZE=31


def load_orig_sizes_and_labels(original_graph_filename):
    "Load the sizes & labels for the original DBG node IDs"
    orig_to_labels = {}
    orig_sizes = {}
    def parse_source_graph_labels(node_id, size, names, vals):
        assert names[0] == 'labels'
        labels = vals[0]
        orig_sizes[node_id] = size
        if labels:
            orig_to_labels[node_id] = list(map(int, labels.strip().split(' ')))

    def nop(*x):
        pass

    with open(original_graph_filename) as fp:
        graph_parser.parse(fp, parse_source_graph_labels, nop)

    return orig_sizes, orig_to_labels


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
    p.add_argument('-q', '--quiet', action='store_true')
    args = p.parse_args()

    basename = os.path.basename(args.catlas_dir)

    original_graph = '%s.gxt' % (basename)
    original_graph = os.path.join(args.catlas_dir, original_graph)
    orig_sizes, orig_to_labels = load_orig_sizes_and_labels(original_graph)

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

    if not all_labels:
        print('** error, no labels', file=sys.stderr)
        sys.exit(1)

    if not args.quiet:
        print('search labels:', search_labels)
        print('all labels:', all_labels)

    found = set()
    found_size = 0
    total_size = 0
    n = 0
    for node_id, subject_mh in mxt_dict.items():
        n += 1
        if n % 100 == 0:
            if not args.quiet:
                print('... search', n, len(mxt_dict))

        size = orig_sizes[node_id]

        total_size += size
        if subject_mh.count_common(query_mh) > 0:
            found.add(node_id)
            found_size += size

    if not args.quiet:
        print('found:', len(found))
        print('found size:', found_size, total_size, round(found_size / total_size, 3))

    pos_nodes = all_nodes.intersection(found)
    neg_nodes = all_nodes - pos_nodes

    def has_search_label(node_id):
        node_labels = orig_to_labels.get(node_id, set())
        return bool(search_labels.intersection(node_labels))

    # true positives: how many nodes did we find that had the right label?
    tp = 0

    # false positives: how many nodes did we find that didn't have right label?
    fp = 0

    for node_id in pos_nodes:
        if has_search_label(node_id):
            tp += orig_sizes[node_id]
        else:
            fp += orig_sizes[node_id]

    # true negatives: how many nodes did we miss that didn't have right label?
    tn = 0

    # false negatives: how many nodes did we miss that did have right label?
    fn = 0

    for node_id in neg_nodes:
        if not has_search_label(node_id):
            tn += orig_sizes[node_id]
        else:
            fn += orig_sizes[node_id]

    if not args.quiet:
        print('')
        print('tp:', tp)
        print('fp:', fp)
        print('fn:', fn)
        print('tn:', tn)
        print('')
    
        try:
            print('sensitivity: %.1f' % (100.0 * tp / (tp + fn))) 
            #If there is no true-positives or false-negatives divides by zero

        except ZeroDivisionError:
            print('Cannot calculate sensitivity')

        try:
            print('specificity: %.1f' % (100.0 * tn / (tn + fp)))
            #If there is no false-positives or true-negatives divides by zero
        except ZeroDivisionError:
            print('Cannot calculate specificity')




if __name__ == '__main__':
    main()
