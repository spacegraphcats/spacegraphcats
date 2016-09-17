#! /usr/bin/env python3

"""
Color a domgraph based on nodes under one or more catlas node.

For example, ::

   color-domgraph.py <catlas> 5 1,3,9 colored.gml

will set the 'colorme' attribute on all domgraph nodes that are underneath
catlas nodes 1, 3, and 9 in the r-5 catlas.
"""
import argparse
import spacegraphcats.graph_parser as parser


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix')
    p.add_argument('catlas_r', type=int)
    p.add_argument('domgraph_nodes')
    p.add_argument('gml_outfile')
    args = p.parse_args()

    assignment_vxt = '%s.assignment.%d.vxt' % (args.catlas_prefix,
                                               args.catlas_r)
    catlas_gxt = '%s.catlas.%d.gxt' % (args.catlas_prefix, args.catlas_r)
    catlas_mxt = '%s.catlas.%d.mxt' % (args.catlas_prefix, args.catlas_r)
    catlas_domgraph = '%s.domgraph.%d.gxt' % (args.catlas_prefix,
                                              args.catlas_r)

    minhash_list = ['tr-1.fa.dump.txt', 'tr-2.fa.dump.txt',
                    'tr-cross.fa.dump.txt']

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

    ### next, find the relevant catlas nodes

    ### next, recurse to find the leaf nodes
    
    def recurse_from(node_id):
        x = []

        beneath = edges.get(node_id, [])
        if beneath:
            for y in beneath:
                x += recurse_from(y)
            return x
        else:
            assert node_id in leaves
            return [leaves_to_domnode[node_id]]

    end_node_to_idx = {}
    domgraph_nodes = args.domgraph_nodes.split(',')
    for n in domgraph_nodes:
        n = int(n)
        for k in recurse_from(n):
            end_node_to_idx[k] = n
    
    print('coloring', len(end_node_to_idx), 'of', len(leaves))

    ### backtrack the leaf nodes to the domgraph
    
    dom_to_orig = {}
    for line in open(assignment_vxt, 'rt'):
        dom_node, orig_list = line.strip().split(',')
        orig_list = list(map(int, orig_list.split(' ')))

        dom_to_orig[int(dom_node)] = orig_list

    x = {}
    for k, idx in end_node_to_idx.items():
        for v in dom_to_orig[k]:
            x[v] = idx
    end_node_to_idx = x

    ### finally, parse the domgraph gxt and convert to a gml w/appropriate attr
    
    with open(args.gml_outfile, 'wt') as fp:
        writer = parser.GmlWriter(fp, directed=False)

        def write_edge(u, v, attributes, values):
            writer.add_edge(u, v, values, attributes)

        def write_vertex(node_id, size, names, vals):
            if node_id in end_node_to_idx:
                names.append('colorme')
                vals.append(str(end_node_to_idx[node_id]))
            writer.add_vertex(node_id, size, vals, names)

        parser.parse(open(catlas_domgraph), write_vertex, write_edge)
    
if __name__ == '__main__':
    main()
