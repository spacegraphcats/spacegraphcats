lint:
	mypy spacegraphcats/*.py search/*.py --ignore-missing-imports

test:
	py.test spacegraphcats search


## Targets:
##
##   search: execute a small build-and-search on the 'acido' data set.
##   15genome-search: execute a medium build-and-search on 15 genomes.

acido-clean:
	-rm -r acido

# build cDBG
acido/cdbg.gxt: data/acido.fa.gz
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 1e9 data/acido.fa.gz

# build catlas
acido/catlas.csv: acido/cdbg.gxt
	python -m spacegraphcats.catlas acido 1

# build minhashes database
acido/minhashes.db: acido/catlas.csv
	python -m search.make_catlas_minhashes acido -k 31 --scaled=1000 --sbt --sigs

# build a search signature
acido/acido-chunk1.fa.gz.sig: data/acido-chunk1.fa.gz
	sourmash compute -k 31 data/acido-chunk1.fa.gz --scaled 500 -f -o acido/acido-chunk1.fa.gz.sig

acido-search: acido/minhashes.db acido/acido-chunk1.fa.gz.sig
	python -m search.search_catlas_with_minhash acido/acido-chunk1.fa.gz.sig acido

acido-frontier-search: acido/minhashes.db acido/acido-chunk1.fa.gz.sig
	python -m search.frontier_search acido/acido-chunk1.fa.gz.sig acido 0.1

acido-frontier-search-optimized: acido/minhashes.db acido/acido-chunk1.fa.gz.sig
	python -m search.frontier_search acido/acido-chunk1.fa.gz.sig acido 0.1  --purgatory


### 

15genome-clean:
	-rm -r 15genome/

# build cDBG
15genome/cdbg.gxt:
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 4e9 data/15genome.fa.gz -o 15genome

# build catlas
15genome/catlas.csv: 15genome/cdbg.gxt
	python -m spacegraphcats.catlas 15genome 3

# build minhashes
15genome/minhashes.db: 15genome/catlas.csv
	python -m search.make_catlas_minhashes -k 31 --scaled=5000 15genome

# run search!
15genome-search: 15genome/minhashes.db
	python -m search.search_catlas_with_minhash data/15genome.5.fa.sig 15genome

15genome-frontier-search: 15genome/minhashes.db
	python -m search.frontier_search data/15genome.5.fa.sig 15genome 0.1

15genome-frontier-search-optimized: 15genome/minhashes.db
	python -m search.frontier_search data/15genome.5.fa.sig 15genome 0.1 --purgatory

####

#
# shewanella.mappedreads.fa is a collection of reads from podar data
# that maps to the Shewanella OS228 genome via bwa aln.  "Real" data,
# with known answer.
#

# prepare reads -- # this is here only for record keeping - never
# needs to be done again.
XXXshewanella.abundtrim.gz:
	trim-low-abund.py --normalize 12 -V -Z 10 -M 2e9 -C 3 -k 21 shewanella.mappedreads.fa -o shewanella.abundtrim.gz --gzip

# download the prepared reads (27 MB) from OSF
shewanella.abundtrim.gz:
	curl -L 'https://osf.io/7az9p/?action=download' > shewanella.abundtrim.gz

# build cDBG
shew-reads/cdbg.gxt: shewanella.abundtrim.gz
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 4e9 shewanella.abundtrim.gz -o shew-reads

# build catlas
shew-reads/catlas.csv: shew-reads/cdbg.gxt
	python -m spacegraphcats.catlas shew-reads 1

# build minhashes
shew-reads/minhashes.db: shew-reads/catlas.csv shew-reads/contigs.fa.gz
	python -m search.make_catlas_minhashes -k 31 --scaled=1000 shew-reads

# download the shewanella genome from OSF
shew-reads/shewanella.fa.gz:
	mkdir -p shew-reads
	curl -L 'https://osf.io/fx4ew/?action=download' > shew-reads/shewanella.fa.gz

# compute shewanella genome signature
shew-reads/shewanella.fa.gz.sig: shew-reads/shewanella.fa.gz
	sourmash compute -k 31 --scaled=1000 shew-reads/shewanella.fa.gz -o shew-reads/shewanella.fa.gz.sig

# run frontier search
shew-search: shew-reads/shewanella.fa.gz.sig shew-reads/minhashes.db
	python -m search.frontier_search shew-reads/shewanella.fa.gz.sig shew-reads 0.1 --purgatory
