#! /usr/bin/env python3
import os
import sys
import argparse
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


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('catlas_r', help='catlas radius', type=int)
    p.add_argument('-x', '--bf-size', type=float, default=1e5)
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

    # load MXT into dict
    graph_minhashes = import_graph_mxt(graphmxt)
    print('imported {} graph minhashes'.format(len(graph_minhashes)))

    factory = GraphFactory(1, args.bf_size, 4)
    tree = SBT(factory)

    print('building tree...')
    for node_id, mh in graph_minhashes.items():
        e = sourmash_lib.Estimators(ksize=MINHASH_K, n=len(mh))
        e.mh = mh
        ss = signature.SourmashSignature('', e, name='{}'.format(node_id))

        leaf = SigLeaf(ss.md5sum(), ss)
        tree.add_node(leaf)

    print('...done with {} minhashes. saving!'.format(len(graph_minhashes)))
    tree.save(args.catlas_prefix)

    sys.exit(0)

if __name__ == '__main__':
    main()
