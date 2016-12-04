#! /usr/bin/env python3
import os
import sys
import argparse
from spacegraphcats import graph_parser
from sourmash_lib import MinHash
from collections import defaultdict
from spacegraphcats.catlas import CAtlas
import time
import sourmash_lib
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import search_minhashes, SigLeaf
from sourmash_lib import signature


MINHASH_SIZE=1000
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
            mh = sourmash_lib.MinHash(len(mins), MINHASH_K)
            for k in mins:
                mh.add_hash(k)
            d[g_id] = mh

    return d


def load_orig_sizes_and_labels(original_graph_filename):
    "Load the sizes & labels for the original DBG node IDs"
    orig_to_labels = defaultdict(list)
    orig_sizes = defaultdict(int)
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
    query_mh = MinHash(len(hashes), MINHASH_K)

    for h in hashes:
        query_mh.add_hash(int(h))
    return query_mh


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('catlas_r', help='catlas radius', type=int)
    p.add_argument('mh_file',
                   help='file containing dumped MinHash signature')
    p.add_argument('--threshold', type=int, default=1)
    p.add_argument('-q', '--quiet', action='store_true')
    p.add_argument('--append-csv', type=str,
                   help='append results in CSV to this file')
    args = p.parse_args()
    
    radius = args.catlas_r
    basename = os.path.basename(args.catlas_prefix)
    graphmxt = '%s.mxt' % (basename,)
    graphmxt = os.path.join(args.catlas_prefix, graphmxt)
    
    catgxt = '%s.catlas.%d.gxt' % (basename, radius)
    catmxt = '%s.catlas.%d.mxt' % (basename, radius)
    catgxt = os.path.join(args.catlas_prefix, catgxt)
    catmxt = os.path.join(args.catlas_prefix, catmxt)

    start = time.time()
    _catlas = CAtlas.read(catgxt, None, args.catlas_r)
    print('loaded CAtlas in {0:.1f} seconds.'.format(time.time() - start))

    dbg_graph = '%s.gxt' % (basename)
    dbg_graph = os.path.join(args.catlas_prefix, dbg_graph)

    orig_sizes, orig_to_labels = \
      load_orig_sizes_and_labels(dbg_graph)

    ### backtrack the leaf nodes to the domgraph

    # here, 'dom_to_orig' is a dictionary mapping domination nodes to
    # a list of original vertices in the De Bruijn graph.
    #
    #   dom_to_orig[dom_node_id] => list of [orig_node_ids]

    assignment_vxt = '%s.assignment.%d.vxt' % (basename, radius)
    assignment_vxt = os.path.join(args.catlas_prefix, assignment_vxt)
    dom_to_orig = load_dom_to_orig(assignment_vxt)

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
    labels_file = os.path.basename(args.catlas_prefix) + '.labels.txt'
    labels_file = os.path.join(args.catlas_prefix, labels_file)
    with open(labels_file, 'rt') as fp:
        labels_names = [ x.strip() for x in fp.readlines() ]

    label_list = list(range(1, len(labels_names) + 1))

    sbt_name = os.path.basename(args.catlas_prefix) + '.' + str(radius)
    print('loading sbt "{}"'.format(sbt_name))
    tree = SBT.load(sbt_name, leaf_loader=SigLeaf.load)

    def search(node, mh, count):
        mins = mh.get_mins()

        if isinstance(node, SigLeaf):
            matches = node.data.estimator.mh.count_common(mh)
        else:  # Node or Leaf, Nodegraph by minhash comparison
            matches = sum(1 for value in mins if node.data.get(value))

        if matches >= count:
            return 1
        return 0

    print('starting searches!')

    siglist = sourmash_lib.signature.load_signatures(args.mh_file)
    sigdict = dict( [ (sig.name(), sig) for sig in siglist ] )

    for label in label_list:
        label_seq_name = labels_names[label - 1]
        print('searching for', label_seq_name, ' ', label)
        sig = sigdict[label_seq_name]

        query_mh = sig.estimator.mh

        leaves = set()
        mins = []
        for leaf in tree.find(search, query_mh, args.threshold):
            node_id = int(leaf.data.name())
            mins.extend(leaf.data.estimator.mh.get_mins())
            leaves.add(node_id)

        from sourmash_lib import MinHash
        mh = MinHash(len(query_mh), 31)
        for k in mins:
            mh.add_hash(k)

        print('XXX', mh.compare(query_mh), query_mh.compare(mh))

        ### count the matches/mismatches between MinHash-found nodes
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

        for node_id in neg_nodes:
            if not has_search_label(node_id):
                tn += dom_sizes[node_id]
            else:
                fn += dom_sizes[node_id]

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
        print('%s - sensitivity: %.1f / specificity / %.1f' % (label_seq_name, sens, spec))

        ## some double checks.

        # all/pos/neg nodes at the end are dominating set members.
        assert not pos_nodes - set(dom_sizes.keys())
        assert not neg_nodes - set(dom_sizes.keys())
        assert not all_nodes - set(dom_sizes.keys())

        if args.append_csv:
            write_header = False
            if not os.path.exists(args.append_csv):
                write_header = True
            with open(args.append_csv, 'at') as outfp:
                if write_header:
                    outfp.write('sens, spec, tp, fp, fn, tn, strategy, searchlevel\n')
                outfp.write('%.1f, %.1f, %d, %d, %d, %d, %s, %d\n' %\
                         (sens, spec, tp, fp, fn, tn, 'sbt',
                          0))

    sys.exit(0)

if __name__ == '__main__':
    main()
