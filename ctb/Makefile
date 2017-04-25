all: acido acido-sbt

clean:
	rm -fr acido

acido:
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 1e9 ../data/acido.fa.gz --scaled 1000 --assemble
	python -m spacegraphcats.catlas acido/acido.gxt acido/acido.catlas acido/acido.domfile 1

acido-sbt:
	python make-sbt-of-mxt.py acido

acido.fa.gz.sig:
	~/dev/sourmash/sourmash compute -k 31 ../data/acido.fa.gz --scaled 100 -f

search: acido acido-sbt acido.fa.gz.sig
	~/dev/sourmash/sourmash sbt_gather acido.sbt.json acido.fa.gz.sig --threshold-bp=100

15genome:
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 4e9 ../data/15genome.fa.gz -o 15genome --scaled 1000 --assemble
	python -m spacegraphcats.catlas 15genome/15genome.gxt 15genome/15genome.catlas 15genome/15genome.domfile 3

