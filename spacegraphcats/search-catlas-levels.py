#! /usr/bin/env python
import argparse
import os
from khmer import MinHash

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('catlas_dir')
    parser.add_argument('-r', '--catlas-radius', type=int, default=5)
    parser.add_argument('signature_txt_file')
    args = parser.parse_args()

    mxtfile = '.catlas.%d.mxt' % args.catlas_radius
    mxtfile = os.path.basename(args.catlas_dir) + mxtfile
    mxtfile = os.path.join(args.catlas_dir, mxtfile)
    
    gxtfile = '.catlas.%d.gxt' % args.catlas_radius
    gxtfile = os.path.basename(args.catlas_dir) + gxtfile
    gxtfile = os.path.join(args.catlas_dir, gxtfile)

    nodes_by_level = {}
    print('reading gxtfile', gxtfile)
    gxt_nodes = [ x.strip().split(',') for x in open(gxtfile) ]
    gxt_nodes = [ x for x in gxt_nodes[1:] if len(x) == 4 ]
    gxt_nodes = [ (int(x[0]), int(x[3])) for x in gxt_nodes ]
    levels = set([ level for (node, level) in gxt_nodes ])
    for level in levels:
        nodes_by_level[level] = set()

    for node, level in gxt_nodes:
        nodes_by_level[level].add(node)

    print('reading mxt file', mxtfile)
        
    mxt_dict = {}
    for line in open(mxtfile):
        node, hashes = line.strip().split(',')
        node = int(node)
        hashes = [ int(h) for h in hashes.split(' ') ]
        mh = MinHash(len(hashes), 31)
        for h in hashes:
            mh.add_hash(h)
            
        mxt_dict[int(node)] = mh

    print('done!')
    print('loading:', args.signature_txt_file)
    data = open(args.signature_txt_file).read().strip()
    hashes = list(map(int, data.split()))
    mh = MinHash(len(hashes), 31)
    for h in hashes:
        mh.add_hash(h)

    for level in range(0, max(nodes_by_level) + 1):
        results = []
        for node in nodes_by_level[level]:
            results.append((
                mxt_dict[node].compare(mh), mh.compare(mxt_dict[node]), node))
        results.sort()
        score1, score2, node = results[-1]
        print('match @ level', level, '%.3f' % score1, '%.3f' % score2, node)

if __name__ == '__main__':
    main()
