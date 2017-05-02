lint:
	mypy spacegraphcats/*.py --ignore-missing-imports

test:
	py.test spacegraphcats


## Targets:
##
##   search: execute a small build-and-search on the 'acido' data set.
##   15genome-search: execute a medium build-and-search on 15 genomes.

acido-clean:
	-rm -r acido

# build cDBG
acido/acido.gxt: data/acido.fa.gz
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 1e9 data/acido.fa.gz

# build catlas
acido/acido.catlas: acido/acido.gxt
	python -m spacegraphcats.catlas acido 1

# build minhashes
acido/acido.minhashes: acido/acido.catlas
	python -m search.make_catlas_minhashes acido -k 31 --scaled=1000 --sbt --sigs

# build a search signature
acido/acido-chunk1.fa.gz.sig: data/acido-chunk1.fa.gz
	sourmash compute -k 31 data/acido-chunk1.fa.gz --scaled 500 -f -o acido/acido-chunk1.fa.gz.sig

acido-search: acido/acido.minhashes acido/acido-chunk1.fa.gz.sig
	python -m search.search_catlas_with_minhash acido/acido-chunk1.fa.gz.sig acido

### 

15genome-clean:
	-rm -r 15genome/

# build cDBG
15genome/15genome.gxt:
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 4e9 data/15genome.fa.gz -o 15genome

# build catlas
15genome/15genome.catlas: 15genome/15genome.gxt
	python -m spacegraphcats.catlas 15genome 3

# build minhashes
15genome/15genome.minhashes: 15genome/15genome.catlas
	python -m search.make_catlas_minhashes -k 31 --scaled=5000 15genome

# run search!
15genome-search: 15genome/15genome.minhashes
	python -m search.search_catlas_with_minhash data/15genome.5.fa.sig 15genome
