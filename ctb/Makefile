all: search

clean:
	rm -fr acido

acido/acido.gxt:
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 1e9 ../data/acido.fa.gz
	python -m spacegraphcats.catlas acido/acido.gxt acido/acido.catlas acido/acido.domfile 1

acido/acido.minhashes: acido/acido.gxt
	python make-sbt-from-contigs.py acido -k 31 --scaled=1000 --sbt

acido.fa.gz.sig:
	sourmash compute -k 31 ../data/acido.fa.gz --scaled 1000 -f

search: acido/acido.minhashes acido.fa.gz.sig
	./search-catlas-with-minhash.py acido.fa.gz.sig acido

15genome-clean:
	rm -r 15genome/

15genome/15genome.gxt:
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 4e9 ../data/15genome.fa.gz -o 15genome

15genome/15genome.catlas: 15genome/15genome.gxt
	python -m spacegraphcats.catlas 15genome/15genome.gxt 15genome/15genome.catlas 15genome/15genome.domfile 3

15genome-sigs: 15genome/15genome.catlas
	./make-sbt-from-contigs.py -k 31 --scaled=5000 15genome

15genome/15genome.minhashes: 15genome-sigs

15genome-search: 15genome/15genome.minhashes
	./search-catlas-with-minhash.py ../data/15genome.5.fa.sig 15genome

15genome-tetra:
	python make-sbt-from-contigs.py -k 4 --scaled=1 --leaves-only 15genome 

15genome-tetra-old:
	extract-long-sequences.py -l 2000 15genome/15genome.gxt.contigs > 15genome/15genome.long-contigs
	sourmash compute -k 4,5 --scaled=1 15genome/15genome.long-contigs --singleton
