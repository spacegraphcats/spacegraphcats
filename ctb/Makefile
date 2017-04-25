all: acido acido-sbt

clean:
	rm -fr acido

acido:
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 1e9 ../data/acido.fa
	python -m spacegraphcats.catlas acido/acido.gxt acido/acido.catlas acido/acido.domgraph 1

acido-sbt:
	python make-sbt-of-mxt.py acido


acido.fa.gz.sig:
	~/dev/sourmash/sourmash compute -k 31 ../data/acido.fa.gz --scaled 100 -f

search:
	~/dev/sourmash/sourmash sbt_gather acido.sbt.json acido.fa.gz.sig --threshold-bp=100
