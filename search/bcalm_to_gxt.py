#! /usr/bin/env python
import screed
import sys
import collections

def main():

    link_d = collections.defaultdict(set)

    gxtfp = open(sys.argv[2], 'wt')
    contigsfp = open(sys.argv[3], 'wt')

    for n, record in enumerate(screed.open(sys.argv[1])):
        if n % 10000 == 0:
            print('...', n, file=sys.stderr)

        name = record.name
        contig_id = int(name.split()[0])
        links = [ x for x in name.split() if x.startswith('L:') ]
        link_ids = [ x.split(':')[2] for x in links ]
        link_ids = [ int(x) for x in link_ids ]

        link_d[contig_id].update(link_ids)

        contigsfp.write('>{}\n{}\n'.format(contig_id, record.sequence))

    gxtfp.write('{}\n'.format(max(link_d.keys()) + 1))

    for k, v in link_d.items():
        for vv in v:
            gxtfp.write('{} {}\n'.format(k, vv))

if __name__ == '__main__':
    main()
