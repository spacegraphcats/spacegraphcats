lint:
	mypy spacegraphcats/*.py search/*.py --ignore-missing-imports

test:
	py.test spacegraphcats search


## Targets:
##
##   acido-search: execute a small build-and-search on the 'acido' data set.
##   15genome-search: execute a medium build-and-search on 15 genomes.
##   shew-search: execute a small build-and-search on a real read data set
##
##   podar-search: execute a search on the full podar data set
##       (requires 8GB RAM)
##   podar-download: download and set up a prebuilt podar catlas.
##       (see https://osf.io/h79um/?show=revision)

acido-clean:
	-rm -r acido

# build cDBG
acido/cdbg.gxt: data/acido.fa.gz
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 2e9 data/acido.fa.gz

# build catlas
acido/catlas.csv: acido/cdbg.gxt
	python -m spacegraphcats.catlas acido 1

# build minhashes database
acido/minhashes_info.json: acido/catlas.csv
	python -m search.make_catlas_minhashes acido -k 31 --scaled=1000

# build a search signature
acido/acido-chunk1.fa.gz.sig: data/acido-chunk1.fa.gz
	sourmash compute -k 31 data/acido-chunk1.fa.gz --scaled 500 -f -o acido/acido-chunk1.fa.gz.sig

acido-simple-search: acido/minhashes_info.json acido/acido-chunk1.fa.gz.sig
	python -m search.search_catlas_with_minhash acido/acido-chunk1.fa.gz.sig acido

acido-frontier-search: acido/minhashes_info.json acido/acido-chunk1.fa.gz.sig
	python -m search.frontier_search acido/acido-chunk1.fa.gz.sig acido 0.1 --fullstats

acido-frontier-search-optimized: acido/minhashes_info.json acido/acido-chunk1.fa.gz.sig
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
15genome/minhashes_info.json: 15genome/catlas.csv
	python -m search.make_catlas_minhashes -k 31 --scaled=5000 15genome

# run search!
15genome-search: 15genome/minhashes_info.json
	python -m search.search_catlas_with_minhash data/15genome.5.fa.sig 15genome

15genome-frontier-search: 15genome/minhashes_info.json
	python -m search.frontier_search data/15genome.5.fa.sig 15genome 0.1

15genome-frontier-search-optimized: 15genome/minhashes_info.json
	python -m search.frontier_search data/15genome.5.fa.sig 15genome 0.1 --purgatory

####

#
# shewanella.mappedreads.fa is a collection of reads from podar data
# that maps to the Shewanella OS228 genome via bwa aln.  "Real" data,
# with known answer.
#

# prepared reads -- this is here only for record keeping & never
# needs to be done again.
XXXshew-reads.abundtrim.gz:
	trim-low-abund.py --normalize 12 -V -Z 10 -M 2e9 -C 3 -k 21 shewanella.mappedreads.fa -o shew-reads.abundtrim.gz --gzip

# download the prepared reads (27 MB) from OSF
shew-reads.abundtrim.gz:
	curl -L 'https://osf.io/7az9p/?action=download' > shew-reads.abundtrim.gz

# build cDBG
shew-reads/cdbg.gxt: shew-reads.abundtrim.gz
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 4e9 shew-reads.abundtrim.gz -o shew-reads

# build catlas
shew-reads/catlas.csv: shew-reads/cdbg.gxt
	python -m spacegraphcats.catlas shew-reads 1

# build minhashes
shew-reads/minhashes_info.json: shew-reads/catlas.csv shew-reads/contigs.fa.gz
	python -m search.make_catlas_minhashes -k 31 --scaled=1000 shew-reads

# download the shewanella genome from OSF
shew-reads/shewanella.fa.gz:
	mkdir -p shew-reads
	curl -L 'https://osf.io/fx4ew/?action=download' > shew-reads/shewanella.fa.gz

# compute shewanella genome signature
shew-reads/shewanella.fa.gz.sig: shew-reads/shewanella.fa.gz
	sourmash compute -k 31 --scaled=1000 shew-reads/shewanella.fa.gz -o shew-reads/shewanella.fa.gz.sig

# run frontier search
shew-search: shew-reads/shewanella.fa.gz.sig shew-reads/minhashes_info.json
	python -m search.frontier_search shew-reads/shewanella.fa.gz.sig shew-reads 0.1 --purgatory

###

#
# SRR606249.keep.fq.gz is the 'podar' data set in reads - from Shakya et al.,
# 2013.  Here it has been prepared from reads that were first QCed as in
# Awad et al. (unpublished), and then normalized and trimmed like so:
#
#    trim-low-abund.py -k 21 -M 8e9 -C 10 -V --normalize 10
#			SRR606249.pe.qc.fq.gz --gzip -o SRR606249.keep.fq.gz
#

# download the prepared reads - 5.3GB in size.
SRR606249.keep.fq.gz:
	curl -L https://osf.io/45xay/?action=download > SRR606249.keep.fq.gz

# download the prepared catlas/minhashes: 250 MB.
podar-download:
	curl -L https://osf.io/g6n4k/?action=download > podar-2017.05.06b.tar.gz
	tar xzf podar-2017.05.06b.tar.gz
	touch podar.ng SRR606249.keep.fq.gz podar/*

# load reads into a nodegraph (8 GB in size)
podar.ng: SRR606249.keep.fq.gz
	load-graph.py -n -M 8e9 -k 31 podar.ng SRR606249.keep.fq.gz

podar/cdbg.gxt: podar.ng SRR606249.keep.fq.gz
	python -m spacegraphcats.build_contracted_dbg -l podar.ng \
		SRR606249.keep.fq.gz -o podar

podar/catlas.csv: podar/cdbg.gxt
	python -m spacegraphcats.catlas podar 3

podar/minhashes.db: podar/cdbg.gxt podar/catlas.csv
	python -m search.make_catlas_minhashes podar -k 31 --scaled=10000

podar-search: podar/minhashes.db
	time python -m search.frontier_search data/mircea-sigs/mircea-rm18.0.fa.sig podar 0.1 --purgatory

####

# build cDBG
twofoo/cdbg.gxt:
	python -m spacegraphcats.build_contracted_dbg -k 21 -M 4e9 shew-reads.abundtrim.gz akker-reads.abundtrim.gz -o twofoo

# build catlas
twofoo/catlas.csv: twofoo/cdbg.gxt
	python -m spacegraphcats.catlas twofoo 1

# build minhashes
twofoo/minhashes.db: twofoo/catlas.csv twofoo/contigs.fa.gz
	python -m search.make_catlas_minhashes -k 21 --scaled=1000 twofoo

twofoo.labels.sqlite: twofoo/catlas.csv
	python -m search.label_cdbg_sqlite twofoo twofoo.fq.gz.bgz twofoo.labels -k 21 -M 1e9

twofoo-extract-1: twofoo/minhashes.db twofoo.labels
	python -m search.extract_reads_by_frontier_sqlite 63-os223.sig twofoo 0.2 -k 21 twofoo.fq.gz.bgz twofoo.labels twofoo.frontier.63.fq

twofoo-extract-p: twofoo/minhashes.db twofoo.labels
	python -m search.extract_contigs_by_frontier 63-os223.sig twofoo 0.2 -k 21 twofoo.frontier.contigs.63.fq

twofoo-extract-bulk:
	python -m search.frontier_search_batch twofoo twofoo.fq.gz.bgz twofoo.labels 2-akker.sig 47-os185.sig 63-os223.sig  -k 21 --savedir foo -o foo/results.csv

twofoo-extract: twofoo/minhashes.db twofoo.labels
	python -m search.extract_reads_by_frontier 63-os223.sig twofoo 0.2 -k 21 twofoo.fq.gz twofoo.labels twofoo.frontier.63.fq
	python -m search.extract_reads_by_frontier 47-os185.sig twofoo 0.2 -k 21 twofoo.fq.gz twofoo.labels twofoo.frontier.47.fq
	python -m search.extract_reads_by_frontier 2-akker.sig twofoo 0.2 -k 21 twofoo.fq.gz twofoo.labels twofoo.frontier.2.fq

twofoo-extract-sqlite: twofoo/minhashes.db twofoo.labels
	python -m search.extract_reads_by_frontier_sqlite 63-os223.sig twofoo 0.2 -k 21 twofoo.labels twofoo.frontier.sql.63.fq
	python -m search.extract_reads_by_frontier_sqlite 47-os185.sig twofoo 0.2 -k 21 twofoo.labels twofoo.frontier.sql.47.fq
	python -m search.extract_reads_by_frontier_sqlite 2-akker.sig twofoo 0.2 -k 21 twofoo.labels twofoo.frontier.sql.2.fq

twofoo.frontier.63.fq.sig: twofoo.frontier.63.fq
	sourmash compute -k 21,31,51 twofoo.frontier.sql.63.fq --scaled=1000 -f

twofoo.frontier.47.fq.sig: twofoo.frontier.47.fq
	sourmash compute -k 21,31,51 twofoo.frontier.sql.47.fq --scaled=1000 -f

twofoo.frontier.2.fq.sig: twofoo.frontier.2.fq
	sourmash compute -k 21,31,51 twofoo.frontier.sql.2.fq --scaled=1000 -f

twofoo-sigs: twofoo.frontier.63.fq.sig twofoo.frontier.2.fq.sig twofoo.frontier.47.fq.sig

twofoo-search: twofoo/minhashes.db
	python -m search.frontier_search shewshew/63.fa.sig twofoo 0.0 --purgatory -k 31 --fullstats
	python -m search.frontier_search shewshew/63.fa.sig twofoo 0.1 --purgatory -k 31 --fullstats
	python -m search.frontier_search shewshew/63.fa.sig twofoo 0.2 --purgatory -k 31 --fullstats
	python -m search.frontier_search shewshew/63.fa.sig twofoo 0.3 --purgatory -k 31 --fullstats
	python -m search.frontier_search shewshew/63.fa.sig twofoo 0.4 --purgatory -k 31 --fullstats
	python -m search.frontier_search shewshew/63.fa.sig twofoo 0.5 --purgatory -k 31 --fullstats
	python -m search.frontier_search shewshew/63.fa.sig twofoo 0.6 --purgatory -k 31 --fullstats
	python -m search.frontier_search shewshew/63.fa.sig twofoo 0.7 --purgatory -k 31 --fullstats
	python -m search.frontier_search shewshew/63.fa.sig twofoo 0.8 --purgatory -k 31 --fullstats
	python -m search.frontier_search shewshew/63.fa.sig twofoo 0.9 --purgatory -k 31 --fullstats
	python -m search.frontier_search shewshew/63.fa.sig twofoo 1.0 --purgatory -k 31 --fullstats
	python -m search.frontier_search shewshew/63.fa.sig twofoo 2.0 --purgatory -k 31 --fullstats

# build cDBG
akker-reads/cdbg.gxt:
	python -m spacegraphcats.build_contracted_dbg -k 31 -M 4e9 akker-reads.abundtrim.gz -o akker-reads

# build catlas
akker-reads/catlas.csv: akker-reads/cdbg.gxt
	python -m spacegraphcats.catlas akker-reads 1

# build minhashes
akker-reads/minhashes.db: akker-reads/catlas.csv akker-reads/contigs.fa.gz
	python -m search.make_catlas_minhashes -k 21 --scaled=1000 akker-reads

akker-reads.abundtrim.gz.bgz: akker-reads.abundtrim.gz
	python -m search.make_bgzf akker-reads.abundtrim.gz

# build reverse index into reads
akker-reads.labels.sqlite: akker-reads/catlas.csv akker-reads/contigs.fa.gz \
		akker-reads.abundtrim.gz.bgz
	python -m search.label_cdbg_sqlite akker-reads \
			akker-reads.abundtrim.gz.bgz akker-reads.labels.sqlite -k 21

# build reverse index into reads
akker-reads.frontier: akker-reads/minhashes.db
	python -m search.frontier_search 2-akker.sig akker-reads 0.2 \
		--fullstats -k 21 --purgatory

# build reverse index into reads
akker-reads.frontier.2.fq: akker-reads.labels.sqlite akker-reads/minhashes.db
	python -m search.extract_reads_by_frontier_sqlite 2-akker.sig akker-reads \
		0.2 -k 21 akker-reads.abundtrim.gz.bgz akker-reads.labels.sqlite \
		akker-reads.frontier.2.fq
