#! /usr/bin/env python
import argparse
import os
from khmer import MinHash

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('catlas_dir')
    parser.add_argument('signature_txt_files', nargs='+')
    parser.add_argument('-l', '--level', type=int, default=None)
    args = parser.parse_args()

    mxtfile = os.path.basename(args.catlas_dir) + '.catlas.5.mxt'
    mxtfile = os.path.join(args.catlas_dir, mxtfile)
    
    gxtfile = os.path.basename(args.catlas_dir) + '.catlas.5.gxt'
    gxtfile = os.path.join(args.catlas_dir, gxtfile)

    gxt_nodes = set()
    if args.level is not None:
        print('reading gxtfile', gxtfile, 'for level', args.level)
        gxt_nodes = [ x.strip().split(',') for x in open(gxtfile) ]
        gxt_nodes = [ x for x in gxt_nodes[1:] if len(x) == 4 ]
        gxt_nodes = [ (int(x[0]), int(x[3])) for x in gxt_nodes ]
        gxt_nodes = set([ a for (a,b) in gxt_nodes if b == args.level ])

        print('read %d nodes' % len(gxt_nodes))

    print('reading mxt file')
        
    mxt_dict = {}
    for line in open(mxtfile):
        node, hashes = line.strip().split(',')
        node = int(node)
        if gxt_nodes and node not in gxt_nodes:
            continue
        hashes = [ int(h) for h in hashes.split(' ') ]
        mh = MinHash(len(hashes), 31)
        for h in hashes:
            mh.add_hash(h)
            
        mxt_dict[int(node)] = mh

    for filename in args.signature_txt_files:
        print('loading:', filename)
        data = open(filename).read().strip()
        hashes = list(map(int, data.split()))
        mh = MinHash(len(hashes), 31)
        for h in hashes:
            mh.add_hash(h)

        print('comparing!')
        print('')
        results = []
        for v in mxt_dict:
            results.append((
                mxt_dict[v].compare(mh), mh.compare(mxt_dict[v]), v))
        results.sort()
        print('sig:', filename)
        print('top N matches at level %d --' % args.level)
        for score1, score2, node in results[-4:]:
            if score1 > 0.000:
                print('%.3f' % score1, '%.3f' % score2, node)

        print('sum: %.3f' % (sum( [ x[0] for x in results ] )))
        print('---')

if __name__ == '__main__':
    main()
