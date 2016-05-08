#! /bin/bash
set -x
set -e
./spacegraphcats/walk-dbg.py -x 5e7 data/tr-cross.fa -o tr-cross.gxt
./spacegraphcats/walk-dbg.py -x 5e7 data/tr-cross.fa -o tr-cross.gml --gml
./spacegraphcats/walk-dbg.py -x 1e8 data/acido-short.fa -o acido-short.gxt
