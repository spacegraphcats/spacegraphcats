#! /usr/bin/env python3
"""
"""
from __future__ import print_function
import argparse
from spacegraphcats import graph_parser
from spacegraphcats.catlas import CAtlas
from spacegraphcats.minhash import MinHash
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

    with open(original_graph_filename) as fp:
        graph_parser.parse(fp, parse_source_graph_labels)

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


def load_mh_dump(mh_filename):
    "Load a MinHash from a file created with 'sourmash dump'."
    with open(mh_filename) as fp:
        hashes = fp.read().strip().split(' ')
    query_mh = MinHash(len(hashes))

    for h in hashes:
        query_mh.add(int(h))
    return query_mh


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('catlas_r', type=int, help='catlas radius to load')
    p.add_argument('mh_file', help='file containing dumped MinHash signature')
    p.add_argument('label_list', type=str,
                   help='list of labels that should correspond to MinHash')
    p.add_argument('-q', '--quiet', action='store_true')
    args = p.parse_args()

    radius = args.catlas_r

    ### first, parse the catlas gxt
    basename = os.path.basename(args.catlas_prefix)
    catgxt = '%s.catlas.%d.gxt' % (basename, radius)
    catmxt = '%s.catlas.%d.mxt' % (basename, radius)
    catgxt = os.path.join(args.catlas_prefix, catgxt)
    catmxt = os.path.join(args.catlas_prefix, catmxt)
    catlas = CAtlas.read(catgxt, catmxt, args.catlas_r)

    print("Loaded catlas with {} leaves".format(len(catlas.leaves())))

    ### get the labels from the original graph

    # here, 'orig_to_labels' is a dictionary mapping De Bruijn graph node IDs
    # to a list of labels for each node.
    #
    #   orig_to_labels[dbg_node_id] => list of [label_ids]

    graphgxt = '%s.gxt' % (basename)
    graphgxt = os.path.join(args.catlas_prefix, graphgxt)
    orig_to_labels = load_orig_to_labels(graphgxt)

    ### backtrack the leaf nodes to the domgraph

    # here, 'dom_to_orig' is a dictionary mapping domination nodes to
    # a list of original vertices in the De Bruijn graph.
    #
    #   dom_to_orig[dom_node_id] => list of [orig_node_ids]

    assignmentvxt = '%s.assignment.%d.vxt' % (basename, radius)
    assignmentvxt = os.path.join(args.catlas_prefix, assignmentvxt)
    dom_to_orig = load_dom_to_orig(assignmentvxt)

    ## now, transfer labels to dom nodes

    dom_labels = {}
    for k, vv in dom_to_orig.items():
        x = dom_labels.get(k, set())
        for v in vv:
            x.update(orig_to_labels.get(v, []))
        dom_labels[k] = x

    ### load mh

    if not args.quiet:
        print('reading mh file', args.mh_file)

    query_mh = load_mh_dump(args.mh_file)

    ### serach the CAtlas!

    leaves = catlas.query(query_mh, 0, CAtlas.Scoring.height_jaccard, CAtlas.Selection.highest_smallest_intersection, CAtlas.Refinement.greedy_coverage)

    if not args.quiet:
        print('found %d leaves as query result' % (len(leaves)))

    ### finally, count the matches/mismatches between MinHash-found nodes
    ### and expected labels.

    search_labels = set([ int(i) for i in args.label_list.split(',') ])
    all_labels = set()
    for vv in dom_labels.values():
        all_labels.update(vv)


    if not args.quiet:
        print("labels we're looking for:", search_labels)
        print("all labels:", all_labels)

    # all_nodes is set of all labeled node_ids on original graph:
    all_nodes = set(dom_labels.keys())

    # pos_nodes is set of MH-matching node_ids from dominating set.
    pos_nodes = set(leaves)

    # neg_nodes is set of non-MH-matching node IDs
    neg_nodes = all_nodes - pos_nodes


    if not args.quiet:
        print('')
        print('pos dom nodes:', len(pos_nodes))
        print('neg dom nodes:', len(neg_nodes))
        print('all dom nodes:', len(all_nodes))

    def has_search_label(node_id):
        node_labels = dom_labels.get(node_id, set())
        return bool(search_labels.intersection(node_labels))

    # true positives: how many nodes did we find that had the right label?
    tp = 0

    # false positives: how many nodes did we find that didn't have right label?
    fp = 0
    
    for dom_node_id in pos_nodes:
        if has_search_label(dom_node_id):
            tp += 1
        else:
            fp += 1

    # true negatives: how many nodes did we miss that didn't have right label?
    tn = 0
    
    # false negatives: how many nodes did we miss that did have right label?
    fn = 0
    
    for dom_node_id in neg_nodes:
        if not has_search_label(dom_node_id):
            tn += 1
        else:
            fn += 1

    if not args.quiet:
        print('')
        print('tp:', tp)
        print('fp:', fp)
        print('fn:', fn)
        print('tn:', tn)
        print('')
    print('sensitivity: %.1f' % (100.0 * tp / (tp + fn)))
    print('specificity: %.1f' % (100.0 * tn / (tn + fp)))

    assert tp + fp + fn + tn == len(all_nodes)

    ## some double checks.

    # all/pos/neg nodes at the end are dominating set members.
    assert not pos_nodes - set(dom_to_orig.keys())
    assert not neg_nodes - set(dom_to_orig.keys())
    assert not all_nodes - set(dom_to_orig.keys())

    sys.exit(0)


if __name__ == '__main__':
    main()
