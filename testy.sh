set -e

mkdir -p ttt
spacegraphcats/walk-dbg.py -x 1e7 data/tr-cross.fa -o ttt/ttt.gxt
spacegraphcats/build-catlas.py ttt 1
spacegraphcats/build-catlas.py ttt 2
spacegraphcats/build-catlas.py ttt 3

spacegraphcats/gxt_to_gml.py ttt/ttt.domgraph.1.gxt ttt/ttt.domgraph.1.gml
spacegraphcats/gxt_to_gml.py ttt/ttt.domgraph.2.gxt ttt/ttt.domgraph.2.gml
spacegraphcats/gxt_to_gml.py ttt/ttt.domgraph.3.gxt ttt/ttt.domgraph.3.gml

spacegraphcats/compute-mh-txt.py data/tr-1.fa -o ttt/tr-1.fa.dump.txt
spacegraphcats/compute-mh-txt.py data/tr-2.fa -o ttt/tr-2.fa.dump.txt
spacegraphcats/compute-mh-txt.py data/tr-cross.fa -o ttt/tr-cross.fa.dump.txt

spacegraphcats/search-catlas-levels.py -r 1 ttt ttt/tr-1.fa.dump.txt
spacegraphcats/search-catlas-levels.py -r 2 ttt ttt/tr-1.fa.dump.txt
spacegraphcats/search-catlas-levels.py -r 3 ttt ttt/tr-1.fa.dump.txt

spacegraphcats/search-catlas-levels.py -r 1 ttt ttt/tr-2.fa.dump.txt
spacegraphcats/search-catlas-levels.py -r 1 ttt ttt/tr-cross.fa.dump.txt

spacegraphcats/search-catlas-mh.py -r 3 ttt ttt/*.dump.txt -l 1
