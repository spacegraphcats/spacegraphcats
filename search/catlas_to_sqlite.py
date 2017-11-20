#! /usr/bin/env python
"""
Convert catlas DAG from CSV to sqlite for searches.
"""
import argparse
import os
import sys

import sqlite3
import gc

from . import search_utils

def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    args = p.parse_args()

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')
    dbfilename = os.path.join(args.catlas_prefix, 'catlas.csv.sqlite')

    if os.path.exists(dbfilename):
        print('removing existing db {}'.format(dbfilename))
        os.unlink(dbfilename)

    print('loading DAG from text files...')
    top_node_id, dag, dag_up, dag_levels = search_utils.load_dag(catlas)

    db = sqlite3.connect(dbfilename)
    cursor = db.cursor()
    cursor.execute('CREATE TABLE top_node (node_id INTEGER)')
    cursor.execute('CREATE TABLE children (node_id INTEGER, child_id INTEGER)');
    cursor.execute('CREATE TABLE first_doms (node_id INTEGER, cdbg_id INTEGER)')
    db.commit()

    cursor.execute('INSERT INTO top_node (node_id) VALUES (?)', (top_node_id,))

    print('...turning it into a SQLite database.')
    for node_id, children in dag.items():
        for child_id in children:
            cursor.execute('INSERT INTO children (node_id, child_id) VALUES (?, ?)', (node_id, child_id))
    db.commit()

    del dag, dag_up, dag_levels
    gc.collect()

    print('loading layer1 from text files...')
    layer1_to_cdbg = search_utils.load_layer1_to_cdbg(catlas, domfile)

    print('...turning it into a sqlite database.')
    for leaf_id, cdbg_list in layer1_to_cdbg.items():
        for cdbg_id in cdbg_list:
            cursor.execute('INSERT INTO first_doms (node_id, cdbg_id) VALUES (?, ?)', (node_id, cdbg_id))
    db.commit()

    sys.exit(0)


if __name__ == '__main__':
    main()
