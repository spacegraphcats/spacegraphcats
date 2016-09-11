#! /usr/bin/env python3
"""
"""
from __future__ import print_function
import argparse
from spacegraphcats import graph_parser
from spacegraphcats.catlas_reader import CAtlasReader
from sourmash_lib import MinHash
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
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('catlas_r', type=int, help='catlas radius to load')
    p.add_argument('mh_files', nargs='+',
                   help='files containing dumped MinHash signature')
    p.add_argument('--strategy', type=str, default='bestnode',
                   help='search strategy: bestnode, xxx')
    p.add_argument('--searchlevel', type=int, default=0)
    p.add_argument('-q', '--quiet', action='store_true')
    p.add_argument('--append-csv', type=str,
                   help='append results in CSV to this file')
    args = p.parse_args()

    ### first, parse the catlas gxt

    catlas = CAtlasReader(args.catlas_prefix, args.catlas_r)

    ### get the labels from the original graph

    # here, 'orig_to_labels' is a dictionary mapping De Bruijn graph node IDs
    # to a list of labels for each node.
    #
    #   orig_to_labels[dbg_node_id] => list of [label_ids]

    orig_sizes, orig_to_labels = \
      load_orig_sizes_and_labels(catlas.original_graph)

    ### backtrack the leaf nodes to the domgraph

    # here, 'dom_to_orig' is a dictionary mapping domination nodes to
    # a list of original vertices in the De Bruijn graph.
    #
    #   dom_to_orig[dom_node_id] => list of [orig_node_ids]

    dom_to_orig = load_dom_to_orig(catlas.assignment_vxt)

    ## now, transfer labels to dom nodes

    dom_labels = {}
    for k, vv in dom_to_orig.items():
        x = dom_labels.get(k, set())
        for v in vv:
            x.update(orig_to_labels.get(v, []))
        dom_labels[k] = x

    ### transfer sizes to dom nodes

    dom_sizes = {}
    for k, vv in dom_to_orig.items():
        x = dom_sizes.get(k, 0)
        for v in vv:
            x += orig_sizes[v]
        dom_sizes[k] = x

    ### load mxt

    # 'mxt_dict' is a dictionary mapping catlas node IDs to MinHash
    # objects.

    mxt_dict = load_mxt_dict(catlas.catlas_mxt)

    # get all labels in graph
    all_labels = set()
    for vv in dom_labels.values():
        all_labels.update(vv)

    label = 1
    for mh_file in args.mh_files:
        ### load search mh
        query_mh = load_mh_dump(mh_file)

        ### next, find the relevant catlas nodes using the MinHash.
        if args.strategy == 'bestnode':
            match_nodes = catlas.find_matching_nodes_best_match(query_mh, mxt_dict)
        elif args.strategy == 'searchlevel':
            match_nodes = catlas.find_matching_nodes_search_level(query_mh, mxt_dict, args.searchlevel)
        elif args.strategy == 'gathermins':
            match_nodes = catlas.find_matching_nodes_gather_mins(query_mh, mxt_dict, args.searchlevel)
        elif args.strategy == 'gathermins2':
            match_nodes = catlas.find_matching_nodes_gather_mins2(query_mh, mxt_dict, args.searchlevel)
        else:
            print('\n*** search strategy not understood:', args.strategy)
            sys.exit(-1)

        leaves = set()
        for match_node in match_nodes:
            leaves.update(catlas.find_level0(match_node))

        ### finally, count the matches/mismatches between MinHash-found nodes
        ### and expected labels.

        search_labels = set([ label ])

        # all_nodes is set of all labeled node_ids on original graph:
        all_nodes = set(dom_labels.keys())

        # pos_nodes is set of MH-matching node_ids from dominating set.
        pos_nodes = set(leaves)

        # neg_nodes is set of non-MH-matching node IDs
        neg_nodes = all_nodes - pos_nodes

        def has_search_label(node_id):
            node_labels = dom_labels.get(node_id, set())
            return bool(search_labels.intersection(node_labels))

        # true positives: how many nodes did we find that had the right label?
        tp = 0

        # false positives: how many nodes did we find that didn't have right label?
        fp = 0

        for dom_node_id in pos_nodes:
            if has_search_label(dom_node_id):
                tp += dom_sizes[dom_node_id]
            else:
                fp += dom_sizes[dom_node_id]

        # true negatives: how many nodes did we miss that didn't have right label?
        tn = 0

        # false negatives: how many nodes did we miss that did have right label?
        fn = 0

        for dom_node_id in neg_nodes:
            if not has_search_label(dom_node_id):
                tn += dom_sizes[dom_node_id]
            else:
                fn += dom_sizes[dom_node_id]

        sens = (100.0 * tp / (tp + fn))
        spec = (100.0 * tn / (tn + fp))
        print('%s - sensitivity: %.1f / specificity / %.1f' % (mh_file, sens, spec))

        ## some double checks.

        # all/pos/neg nodes at the end are dominating set members.
        assert not pos_nodes - set(dom_to_orig.keys())
        assert not neg_nodes - set(dom_to_orig.keys())
        assert not all_nodes - set(dom_to_orig.keys())

        if args.append_csv:
            write_header = False
            if not os.path.exists(args.append_csv):
                write_header = True
            with open(args.append_csv, 'at') as outfp:
                if write_header:
                    outfp.write('sens, spec, tp, fp, fn, tn, strategy, searchlevel\n')
                outfp.write('%.1f, %.1f, %d, %d, %d, %d, %s, %d\n' %\
                         (sens, spec, tp, fp, fn, tn, args.strategy,
                          args.searchlevel))
        label += 1

    sys.exit(0)


if __name__ == '__main__':
    main()
