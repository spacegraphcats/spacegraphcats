#! /usr/bin/env python
import argparse
import os
from khmer import MinHash


KSIZE=31


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('catlas_dir')
    parser.add_argument('-r', '--catlas-radius', type=int, default=5)
    parser.add_argument('signature_txt_files', nargs='+')
    parser.add_argument('-l', '--level', type=int, default=None)
    args = parser.parse_args()

    mxtfile = '.catlas.%d.mxt' % args.catlas_radius
    mxtfile = os.path.basename(args.catlas_dir) + mxtfile
    mxtfile = os.path.join(args.catlas_dir, mxtfile)
    
    gxtfile = '.catlas.%d.gxt' % args.catlas_radius
    gxtfile = os.path.basename(args.catlas_dir) + gxtfile
    gxtfile = os.path.join(args.catlas_dir, gxtfile)

    gxt_nodes = set()
    gxt_levels = {}
    if args.level is not None:
        print('reading gxtfile', gxtfile, 'for level', args.level)
        gxt_nodes = [ x.strip().split(',') for x in open(gxtfile) ]
        gxt_nodes = [ x for x in gxt_nodes[1:] if len(x) == 4 ]
        gxt_nodes = [ (int(x[0]), int(x[3])) for x in gxt_nodes ]
        gxt_nodes = set([ a for (a,b) in gxt_nodes if b == args.level ])

        print('read %d nodes' % len(gxt_nodes))
    else:
        print('searching ALL levels at once')
        gxt_nodes = [ x.strip().split(',') for x in open(gxtfile) ]
        gxt_nodes = [ x for x in gxt_nodes[1:] if len(x) == 4 ]
        gxt_levels = dict([ (int(x[0]), int(x[3])) for x in gxt_nodes ])
        gxt_nodes = set()

    print('reading mxt file')
        
    mxt_dict = {}
    for line in open(mxtfile):
        node, hashes = line.strip().split(',')
        node = int(node)
        if gxt_nodes and node not in gxt_nodes:
            continue
        if args.level and 0:
            fp = open(mxtfile + '.node%d.level%d.mh' % (node, args.level),
                      'wt')
            fp.write(hashes)
            fp.close()
        hashes = [ int(h) for h in hashes.split(' ') ]
        mh = MinHash(len(hashes), KSIZE)
        for h in hashes:
            mh.add_hash(h)
            
        mxt_dict[int(node)] = mh

    for filename in args.signature_txt_files:
        print('loading:', filename)
        data = open(filename).read().strip()
        hashes = list(map(int, data.split()))
        mh = MinHash(len(hashes), KSIZE)
        for h in hashes:
            mh.add_hash(h)

        print('comparing!')
        print('')
        results = []
        for v in mxt_dict:
            level = gxt_levels.get(v, args.level)
            results.append((
                mxt_dict[v].compare(mh), mh.compare(mxt_dict[v]), v, level, mxt_dict[v]))
        results.sort(key=lambda x: x[1])
        print('sig:', filename)
        if args.level:
            print('top N matches at level %d --' % args.level)
        else:
            print('top N matches across all levels')
        for score1, score2, node, level, mxt_mh in results[-4:]:
            if score1 > 0.000:
                print(",".join((filename,'%.3f' % score1, '%.3f' % score2, str(node), str(level))))

        if args.level:
            sum_mh = MinHash(10000000, KSIZE)
            for score1, score2, node, level, mxt_mh in results:
                if score1 > 0.05:
                    sum_mh.merge(mxt_mh)
            print('sum mh: %.3f %.3f' % (sum_mh.compare(mh), mh.compare(sum_mh)))
            #print('sum: %.3f' % (sum( [ x[0] for x in results ] )))
        print('---')

if __name__ == '__main__':
    main()
