#! /usr/bin/env python
import sourmash_lib
import screed
import khmer.utils

both_sig = sourmash_lib.load_one_signature('shew-reads.contigs.sig', select_ksize=31)
reads_sig = sourmash_lib.load_one_signature('shew-reads.abundtrim.gz.sig', select_ksize=31)

diff_mins = set(reads_sig.minhash.get_mins()) - set(both_sig.minhash.get_mins())
#print(diff_mins)

found_mins = set()

outfp = open('shew-reads.nosig.fa', 'wt')
for n, record in enumerate(screed.open('shew-reads.abundtrim.gz')):
    if n % 10000 == 0:
        print('...', n)
    mh = sourmash_lib.MinHash(0, 31, scaled=1000)
    mh.add_sequence(record.sequence, True)
    if diff_mins.intersection(mh.get_mins()):
        found_mins.update(mh.get_mins())
        outfp.write('>{}\n{}\n'.format(record.name, record.sequence))

diff_mins -= found_mins
print('still missing:', diff_mins)
