#! /usr/bin/env python3
import os
import sys
import sqlite3
import argparse
from sourmash_lib import MinHash
from collections import defaultdict
from spacegraphcats.catlas import CAtlas


MINHASH_SIZE=1000
MINHASH_K=0


def connect_db(filepath):
    conn = sqlite3.connect(filepath)
    cur = conn.cursor()
    
    cur.execute("PRAGMA synchronous='OFF'")
    cur.execute("PRAGMA locking_mode=EXCLUSIVE")

    cur.execute('''CREATE TABLE graph_minhashes (node_id INTEGER PRIMARY KEY,
                                                 mins TEXT)''')
    cur.execute('''CREATE TABLE domgraph_minhashes (node_id INTEGER PRIMARY KEY,
                                                 mins TEXT)''')
    cur.execute('''CREATE TABLE catlas_minhashes (node_id INTEGER PRIMARY KEY,
                                                  mins TEXT)''')

    return conn

def import_graph_mxt(conn, mxt):
    n = 0
    with open(mxt, 'rt') as fp:
        cur = conn.cursor()
        
        for line in fp:
            line = line.split(',')
            g_id = int(line[0])
            mins = line[1]

            cur.execute('INSERT INTO graph_minhashes (node_id, mins) VALUES (?, ?)', (g_id, mins))
            n += 1

    conn.commit()

    return n


def export_catlas_mxt(conn, fp):
    n = 0
    cur = conn.cursor()

    cur.execute('SELECT node_id, mins FROM catlas_minhashes')
    for (node_id, mins) in cur:
        fp.write('{},{}\n'.format(node_id, mins))
        n += 1
        
    return n


def merge_nodes(cur, child_node_list, parent_node, from_tablename, to_tablename):
    minlist = []

    for graph_node in child_node_list:
        cur.execute('SELECT mins FROM {} WHERE node_id=?'.format(from_tablename),
                    (graph_node,))
        try:
            mins = cur.fetchone()[0]
        except TypeError as e:
            print('ERROR - no results for node {} in table {} (catlas {})'.format(graph_node, from_tablename, parent_node))
            continue
        mins = list(map(int, mins.split()))
        minlist.append(mins)

    #assert len(minlist) == len(child_node_list)

    # merge (note that minhashes are implicitly merged when you simply add all
    # the hashes to one MH).
    merged_mh = MinHash(MINHASH_SIZE, MINHASH_K)
    for mins in minlist:
        for h in mins:
            merged_mh.add_hash(h)

    # add into catlas minhashes
    merged_mins = " ".join(map(str, merged_mh.get_mins()))
    try:
        cur.execute('INSERT INTO {} (node_id, mins) VALUES (?, ?)'.format(to_tablename),
                    (parent_node, merged_mins))
    except sqlite3.IntegrityError as e:
        print('XXX', child_node_list, catlas_node)
        print(str(e))


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


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('catlas_r', help='catlas radius', type=int)
    args = p.parse_args()
    
    radius = args.catlas_r
    basename = os.path.basename(args.catlas_prefix)
    graphmxt = '%s.mxt' % (basename,)
    graphmxt = os.path.join(args.catlas_prefix, graphmxt)
    
    catgxt = '%s.catlas.%d.gxt' % (basename, radius)
    catmxt = '%s.catlas.%d.mxt' % (basename, radius)
    catgxt = os.path.join(args.catlas_prefix, catgxt)
    catmxt = os.path.join(args.catlas_prefix, catmxt)

    catmxt_db = catmxt + '.db'

    # load mapping between dom nodes and cDBG/graph nodes:
    assignment_vxt = '%s.assignment.%d.vxt' % (basename, radius)
    assignment_vxt = os.path.join(args.catlas_prefix, assignment_vxt)
    dom_to_orig = load_dom_to_orig(assignment_vxt)
    
    # load MXT into sqlite database
    conn = connect_db(catmxt_db)
    n = import_graph_mxt(conn, graphmxt)
    print('imported {} graph minhashes'.format(n))

    # create minhashes for catlas leaf nodes.
    cur = conn.cursor()

    for n, (dom_node, graph_nodes) in enumerate(dom_to_orig.items()):
        merge_nodes(cur, graph_nodes, dom_node, 'graph_minhashes', 'domgraph_minhashes')
    print('created {} leaf node MinHashes via merging'.format(n + 1))

    # this gives us the leaf level minhashes that we need for the rest.
    # now, eliminate the assignment dict & go for the catlas structure.
    del dom_to_orig
    catlas = CAtlas.read(catgxt, None, args.catlas_r)

    # for level 0, merge the shadows (domgraph nodes)
    print('merging level 0')
    n = 0
    m = 0
    select = lambda node: node.level == 0
    for node in catlas.nodes(select):
        domgraph_nodes = node.shadow()
        merge_nodes(cur, domgraph_nodes, node.id,
                    'domgraph_minhashes',
                    'catlas_minhashes')
        n += 1
        m += len(domgraph_nodes)
    print('level 0: merged {} children into {} nodes'.format(n, m))

    # for each level above 0, merge the children
    for level in range(1, catlas.level + 1):
        print('merging at level:', level)
        n = 0
        m = 0
        select = lambda node: node.level == level
        for node in catlas.nodes(select):
            merge_nodes(cur, [ n.id for n in node.children ],
                        node.id, 'catlas_minhashes', 'catlas_minhashes')
            n += 1
            m += len(node.children)
        print('level {}: merged {} children into {} nodes'.format(level, n, m))
    conn.commit()

    with open(catmxt, 'wt') as fp:
        n = export_catlas_mxt(conn, fp)
        print('exported {} catlas minhashes to {}'.format(n, catmxt))

    sys.exit(0)

if __name__ == '__main__':
    main()
