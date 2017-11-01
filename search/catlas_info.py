#! /usr/bin/env python
import argparse
import os
import sys
import leveldb

from .search_utils import (load_dag, load_layer0_to_cdbg, load_minhash)
from .search_utils import get_minhashdb_name


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix')
    p.add_argument('-k', '--ksize', default=31, type=int,
                   help='list of k-mer sizes (default: 31)')
    args = p.parse_args()

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')
    contigfile = os.path.join(args.catlas_prefix, "contigs.fa.gz")

    top_node_id, dag, dag_up, dag_levels = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    print('top catlas node {} has {} children.'.format(top_node_id,
                                                       len(dag[top_node_id])))
    
    layer0_to_cdbg = load_layer0_to_cdbg(catlas, domfile)
    x = set()
    for v in layer0_to_cdbg.values():
        x.update(v)
    print('{} layer 0 catlas nodes, corresponding to {} cDBG nodes.'.format(len(layer0_to_cdbg), len(x)))

    db_path = get_minhashdb_name(args.catlas_prefix, args.ksize, 0, 0)
    if not db_path:
        print('** ERROR, minhash DB does not exist for {}'.format(ksize),
              file=sys.stderr)
        sys.exit(-1)
    minhash_db = leveldb.LevelDB(db_path)

    top_mh = load_minhash(top_node_id, minhash_db)
    print('top node minhash has {} mins, k={} scaled={}'.format(len(top_mh.get_mins()), top_mh.ksize, top_mh.scaled))
    print(' => ~{:g} {}-mers total.'.format(len(top_mh.get_mins()) * top_mh.scaled, top_mh.ksize))

    num_empty = 0
    for child_node in dag[top_node_id]:
        mh = load_minhash(child_node, minhash_db)
        if not mh:
            num_empty += 1

    print('{} empty child nodes (of {}) beneath top node'.format(num_empty, len(dag[top_node_id])))


if __name__ == '__main__':
    # import cProfile
    # cProfile.run('main()', 'search_stats')

    main()
