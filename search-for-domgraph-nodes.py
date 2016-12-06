#! /usr/bin/env python
"""
"""
from __future__ import print_function
import argparse
from spacegraphcats import graph_parser
from spacegraphcats.catlas import CAtlas
from spacegraphcats.graph import VertexDict
from collections import defaultdict
from sourmash_lib import MinHash
import sourmash_lib.signature
import os
import sys
import time
import cProfile
import pstats


KSIZE=31


def load_orig_sizes_and_labels(original_graph_filename):
    "Load the sizes & labels for the original DBG node IDs"
    orig_to_labels = defaultdict(list)
    orig_sizes = defaultdict(int)
    def parse_source_graph_labels(node_id, size, names, vals):
        assert names[0] == 'labels'
        orig_sizes[node_id] = size

        if vals:
            labels = vals[0]
            if labels:
                labels = list(map(int, labels.strip().split(' ')))
                orig_to_labels[node_id] = labels

    def nop(*x):
        pass

    with open(original_graph_filename) as fp:
        graph_parser.parse(fp, parse_source_graph_labels, nop)

    return orig_sizes, orig_to_labels


def load_dom_to_orig(assignment_vxt_file):
    "Load the mapping between level0/dom nodes and the original DBG nodes."
    dom_to_orig = defaultdict(list)
    for line in open(assignment_vxt_file, 'rt'):
        orig_node, dom_list = line.strip().split(',')
        orig_node = int(orig_node)
        dom_list = list(map(int, dom_list.split(' ')))

        for dom_node in dom_list:
            dom_to_orig[dom_node].append(orig_node)

    return dom_to_orig


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
    p.add_argument('mh_file',
                   help='file containing MinHash signatures')
    p.add_argument('--strategy', type=str, default='bestnode',
                   help='search strategy: bestnode, xxx')
    p.add_argument('--searchlevel', type=int, default=0)
    p.add_argument('-q', '--quiet', action='store_true')
    p.add_argument('--append-csv', type=str,
                   help='append results in CSV to this file')
    p.add_argument('--label-list', type=str,
                   help='list of labels that should correspond to MinHash; ' + \
                        'defaults to index of mh_files.')
    args = p.parse_args()

    # load the CAtlas.
    _radius = args.catlas_r
    _basename = os.path.basename(args.catlas_prefix)
    _catgxt = '%s.catlas.%d.gxt' % (_basename, _radius)
    _catmxt = '%s.catlas.%d.mxt' % (_basename, _radius)
    _catgxt = os.path.join(args.catlas_prefix, _catgxt)
    _catmxt = os.path.join(args.catlas_prefix, _catmxt)

    start = time.time()
    _catlas = CAtlas.read(_catgxt, _catmxt, args.catlas_r)
    print('loaded CAtlas in {0:.1f} seconds.'.format(time.time() - start))

    #delim_str = "#"*40
    #print(delim_str)
    #print("Catlas has {} levels".format(_catlas.level))
    #for i,l in enumerate(_catlas.bfs()):
    #    print(i,len(l))
    #print(delim_str)
    
    ### get the labels from the original graph

    # here, 'orig_to_labels' is a dictionary mapping De Bruijn graph node IDs
    # to a list of labels for each node.
    #
    #   orig_to_labels[dbg_node_id] => list of [label_ids]

    _dbg_graph = '%s.gxt' % (_basename)
    _dbg_graph = os.path.join(args.catlas_prefix, _dbg_graph)

    orig_sizes, orig_to_labels = \
      load_orig_sizes_and_labels(_dbg_graph)

    ### backtrack the leaf nodes to the domgraph

    # here, 'dom_to_orig' is a dictionary mapping domination nodes to
    # a list of original vertices in the De Bruijn graph.
    #
    #   dom_to_orig[dom_node_id] => list of [orig_node_ids]

    _assignment_vxt = '%s.assignment.%d.vxt' % (_basename, _radius)
    _assignment_vxt = os.path.join(args.catlas_prefix, _assignment_vxt)
    dom_to_orig = load_dom_to_orig(_assignment_vxt)

    ## now, transfer labels to dom nodes

    dom_labels = defaultdict(set)
    for k, vv in dom_to_orig.items():
        for v in vv:
            dom_labels[k].update(orig_to_labels[v])

    ### transfer sizes to dom nodes

    dom_sizes = defaultdict(int)
    for k, vv in dom_to_orig.items():
        for v in vv:
            dom_sizes[k] += orig_sizes[v]

    end = time.time()
    print('done loading all the things - {0:.1f} seconds.'.format(end - start))

    # get all labels in graph
    all_labels = set()
    for vv in dom_labels.values():
        all_labels.update(vv)

    # construct list of labels to search for
    siglist = sourmash_lib.signature.load_signatures(args.mh_file)
    sigdict = dict( [ (sig.name(), sig) for sig in siglist ] )

    labels_file = os.path.basename(args.catlas_prefix) + '.labels.txt'
    labels_file = os.path.join(args.catlas_prefix, labels_file)
    with open(labels_file, 'rt') as fp:
        labels_names = [ x.strip().split(' ', 1) for x in fp.readlines() ]
        labels_names = dict([ (int(x), sigdict[y]) for (x, y) in labels_names ])

    print('starting searches!')
    for label in labels_names:
        query_sig = labels_names[label]
        print('searching for sequence {} / label {}'.format(query_sig.name(),
                                                            label))

        query_mh = query_sig.estimator.mh

        prof = cProfile.Profile()
        prof.enable()
        q_start = time.time()

        ### next, find the relevant catlas nodes using the MinHash.
        if args.strategy == 'bestnode':
            match_nodes = _catlas.query_best_match(query_mh)
        elif args.strategy == 'searchlevel':
            match_nodes = _catlas.query_level(query_mh, args.searchlevel)
        elif args.strategy == 'gathermins':
            match_nodes = _catlas.query_gather_mins(query_mh, args.searchlevel)
        elif args.strategy == 'gathermins2':
            match_nodes = _catlas.query_gather_mins(query_mh, args.searchlevel, expand=True)
        elif args.strategy == 'frontier-jacc':
            match_nodes  = _catlas.query(query_mh, 0, CAtlas.Scoring.jaccard, CAtlas.Selection.largest_weighted_intersection_cached(), CAtlas.Refinement.greedy_coverage)
        elif args.strategy == 'frontier-jacc-bl':
            match_nodes  = _catlas.query_blacklist(query_mh, 0, CAtlas.Scoring.jaccard, CAtlas.Selection.largest_weighted_intersection(), CAtlas.Refinement.greedy_coverage)
        elif args.strategy == 'frontier-height':
            match_nodes  = _catlas.query(query_mh, 0, CAtlas.Scoring.avg_height, CAtlas.Selection.largest_intersection_height, CAtlas.Refinement.greedy_coverage)
        elif args.strategy == 'frontier-height-bl':
            match_nodes  = _catlas.query_blacklist(query_mh, 0, CAtlas.Scoring.avg_height, CAtlas.Selection.largest_intersection_height, CAtlas.Refinement.greedy_coverage)
        elif args.strategy == 'frontier-worst':
            match_nodes  = _catlas.query(query_mh, 0, CAtlas.Scoring.worst_node, CAtlas.Selection.largest_weighted_intersection, CAtlas.Refinement.greedy_coverage)
        elif args.strategy == 'frontier-worst-bl':
            match_nodes  = _catlas.query_blacklist(query_mh, 0, CAtlas.Scoring.worst_node, CAtlas.Selection.largest_weighted_intersection, CAtlas.Refinement.greedy_coverage)
        else:
            print('\n*** search strategy not understood:', args.strategy)
            sys.exit(-1)

        f_time = time.time()
        print("Frontier searched in {0:.1f} seconds.".format(f_time-q_start))

        leaves = set()
        for n in (match_nodes):
            leaves.update(n.shadow())
        q_end = time.time()

        print("Query completed in {0:.1f} seconds.".format(q_end-q_start))

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
            node_labels = dom_labels[node_id]
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

        if not args.quiet and not args.append_csv:
            print('')
            print('tp:', tp)
            print('fp:', fp)
            print('fn:', fn)
            print('tn:', tn)
            print('')

        sens = 0
        if tp + fn:
            sens = (100.0 * tp / (tp + fn))
        if tn + fp:
            spec = (100.0 * tn / (tn + fp))
        print('%s - sensitivity: %.1f / specificity / %.1f' % \
              (query_sig.name(), sens, spec))

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
        prof.disable()
        ps = pstats.Stats(prof).sort_stats('time')
        ps.print_stats()

    sys.exit(0)


if __name__ == '__main__':
    main()
