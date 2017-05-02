## Targets:
##
##   search: execute a small build-and-search on the 'acido' data set.
##   15genome-search: execute a medium build-and-search on 15 genomes.

all: search

clean:
	-rm -r acido acido.fa.gz.sig

# build cDBG
acido/acido.gxt: ../data/acido.fa.gz
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 1e9 ../data/acido.fa.gz

# build catlas
acido/acido.catlas: acido/acido.gxt
	python -m spacegraphcats.catlas acido/acido.gxt acido/acido.catlas acido/acido.domfile 1

# build minhashes
acido/acido.minhashes: acido/acido.catlas
	python make-catlas-minhashes.py acido -k 31 --scaled=1000 --sbt --sigs

# build a search signature
acido-chunk1.fa.gz.sig: ../data/acido-chunk1.fa.gz
	sourmash compute -k 31 ../data/acido-chunk1.fa.gz --scaled 500 -f

search: acido/acido.minhashes acido-chunk1.fa.gz.sig
	./search-catlas-with-minhash.py acido-chunk1.fa.gz.sig acido

### 

15genome-clean:
	-rm -r 15genome/

# build cDBG
15genome/15genome.gxt:
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 4e9 ../data/15genome.fa.gz -o 15genome

# build catlas
15genome/15genome.catlas: 15genome/15genome.gxt
	python -m spacegraphcats.catlas 15genome/15genome.gxt 15genome/15genome.catlas 15genome/15genome.domfile 3

# build minhashes
15genome/15genome.minhashes: 15genome/15genome.catlas
	./make-catlas-minhashes.py -k 31 --scaled=5000 15genome

# run search!
15genome-search: 15genome/15genome.minhashes
	./search-catlas-with-minhash.py ../data/15genome.5.fa.sig 15genome

15genome-tetra:
	./make-sbt-from-contigs.py -k 4 --scaled=1 --leaves-only 15genome 

15genome-tetra-old:
	extract-long-sequences.py -l 2000 15genome/15genome.gxt.contigs > 15genome/15genome.long-contigs
	sourmash compute -k 4,5 --scaled=1 15genome/15genome.long-contigs --singleton
