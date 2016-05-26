#! /usr/bin/env python
"""
"""
import argparse
import parser                             #spacegraphcats parser
from khmer import MinHash
import os
import sys


KSIZE=31


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix')
    p.add_argument('catlas_r', type=int)
    p.add_argument('mh_file')
    p.add_argument('label_list', type=str)
    args = p.parse_args()

    assignment_vxt = '%s.assignment.%d.vxt' % (args.catlas_prefix,
                                               args.catlas_r)
    assignment_vxt = os.path.join(args.catlas_prefix, assignment_vxt)
    
    catlas_gxt = '%s.catlas.%d.gxt' % (args.catlas_prefix, args.catlas_r)
    catlas_gxt = os.path.join(args.catlas_prefix, catlas_gxt)
    
    catlas_mxt = '%s.catlas.%d.mxt' % (args.catlas_prefix, args.catlas_r)
    catlas_mxt = os.path.join(args.catlas_prefix, catlas_mxt)
                                  
    catlas_domgraph = '%s.domgraph.%d.gxt' % (args.catlas_prefix,
                                              args.catlas_r)
    catlas_domgraph = os.path.join(args.catlas_prefix, catlas_domgraph)

    original_graph = '%s.gxt' % (args.catlas_prefix)
    original_graph = os.path.join(args.catlas_prefix, original_graph)

    ### first, parse the catlas gxt
    
    edges = {}
    vertices = {}
    roots = []
    leaves = []
    leaves_to_domnode = {}

    def add_edge(a, b, *extra):
        x = edges.get(a, [])
        x.append(b)
        edges[a] = x

    def add_vertex(node_id, size, names, vals):
        assert names[0] == 'vertex'
        assert names[1] == 'level'
        vertex, level = vals
        if vertex == 'root':
            roots.append(node_id)

        if level == '0':
            leaves.append(node_id)
            leaves_to_domnode[node_id] = int(vertex)

    parser.parse(open(catlas_gxt), add_vertex, add_edge)

    ### get the labels from the original graph

    orig_to_labels = {}
    def parse_source_graph_labels(node_id, size, names, vals):
        assert names[0] == 'labels'
        labels = vals[0]
        if labels:
            orig_to_labels[node_id] = list(map(int, labels.strip().split(' ')))
        
    def nop(*x):
        pass

    parser.parse(open(original_graph), parse_source_graph_labels, nop)

    ### backtrack the leaf nodes to the domgraph
    
    dom_to_orig = {}
    for line in open(assignment_vxt, 'rt'):
        orig_node, dom_list = line.strip().split(',')
        orig_node = int(orig_node)
        dom_list = list(map(int, dom_list.split(' ')))

        for dom_node in dom_list:
            x = dom_to_orig.get(dom_node, [])
            x.append(orig_node)
            dom_to_orig[dom_node] = x

    ### next, build function to find the leaf nodes
    
    def recurse_from(node_id):
        x = []

        beneath = edges.get(node_id, [])
        if beneath:
            # recurse!
            for y in beneath:
                x += recurse_from(y)
            return x
        else:
            # only thing there should be nothing beneath are the leaves...
            assert node_id in leaves
            # convert leaves into the original domination graph coordinates.
            return [leaves_to_domnode[node_id]]

    ### load mxt

    print('reading mxt file', catlas_mxt)
        
    mxt_dict = {}
    for line in open(catlas_mxt):
        node, hashes = line.strip().split(',')
        node = int(node)
        hashes = [ int(h) for h in hashes.split(' ') ]
        mh = MinHash(len(hashes), KSIZE)
        for h in hashes:
            mh.add_hash(h)
            
        mxt_dict[int(node)] = mh

    ### load mh

    print('reading mh file', args.mh_file)

    hashes = open(args.mh_file).read().strip().split(' ')
    query_mh = MinHash(len(hashes), KSIZE)

    for h in hashes:
        query_mh.add_hash(int(h))
        
    ### next, find the relevant catlas nodes

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
    leaves = set(recurse_from(best_match_node))

    print('found %d domgraph leaves under catlas node %d' % (len(leaves),
                                                      best_match_node))

    ### finally, count!
    
    orig_nodes = []
    for k in leaves:
        orig_nodes += dom_to_orig[k]

    found_label_counts = {}
    for k in orig_nodes:
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
    
    print('looking for labels:', " ".join([str(i) for i in label_list]))
    print('actually found:', found_label_counts)
    print('all label counts:', all_label_counts)

    print('tp:', tp)
    print('fp:', fp)
    print('fn:', fn)
    print('tn:', tn)

    sys.exit(0)


if __name__ == '__main__':
    main()
