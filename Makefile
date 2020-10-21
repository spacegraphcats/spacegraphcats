all: test

flakes:
	# TODO: remove from exceptions F841
	flake8 spacegraphcats/ --ignore E252,E501,E226,E201,E202,266,W504,E266,E127,E302,W293,W291,F841

lint:
	#mypy spacegraphcats/*.py

test:
	pytest spacegraphcats

#
# akker-reads.abundtrim.gz is a collection of reads from podar data
# that maps to the Akkermansia muciniphila ATCC BAA-835 genome via bwa aln.
# "Real" data, with known answer.  There does not appear to be significant
# overlap with other genomes in the Podar data set; so, no significant strain
# variation.
#

akker-reads.abundtrim.gz:
	curl -o akker-reads.abundtrim.gz -L https://osf.io/dk7nb/download

#
# shew-reads.abundtrim.gz is a collection of reads from podar data
# that maps to the Shewanella OS223 genome via bwa aln.  "Real" data,
# with known answer.  Note that there is significant overlap with the
# Shewanella OS185 genome; this is a data set with significant strain
# variation.
#

shew-reads.abundtrim.gz:
	curl -L 'https://osf.io/7az9p/?action=download' > shew-reads.abundtrim.gz

#
# twofoo targets, below, use a synthetic mixture of reads from podar data -
# the shew-reads.abundtrim.gz (mapping to Shewanella baltica OS223) and
# akker-reads.abundtrim.gz (mapping to Akkermansia muciniphila ATCC BAA-835).
# Many of the shew-reads also map to S. baltica OS185, while the akker-reads
# do not; so this is a good mixture for testing the effects of strain variation
# on catlas foo.

# make synthetic mix data set 'twofoo'
twofoo.fq.gz: shew-reads.abundtrim.gz akker-reads.abundtrim.gz
	gunzip -c shew-reads.abundtrim.gz akker-reads.abundtrim.gz | gzip -9c > twofoo.fq.gz

#
# these are targets for quick testing and/or obscure script testing.
#
# dory-test runs through the entire pipeline on a Doryteuthis transcriptome
# subset.
#

dory-test: data/dory-subset.fa data/dory-head.fa
	python -m spacegraphcats dory-test search
	python -m spacegraphcats dory-test extract_reads
	python -m spacegraphcats dory-test extract_contigs
	python -m spacegraphcats dory-test multifasta_query

#twofoo/bcalm.twofoo.k31.unitigs.fa:
#	mkdir -p twofoo
#	curl -L -o twofoo/bcalm.twofoo.k31.unitigs.fa.gz https://osf.io/zp49s/download
#	gunzip twofoo/bcalm.twofoo.k31.unitigs.fa.gz

twofoo-test: twofoo.fq.gz # twofoo/bcalm.twofoo.k31.unitigs.fa
	python -m spacegraphcats twofoo search
	python -m spacegraphcats twofoo hashval_query
	python -m spacegraphcats twofoo extract_reads_for_hashvals
	python -m spacegraphcats.search.characterize_catlas_regions twofoo_k31_r1 twofoo_k31_r1.vec
	python -m spacegraphcats twofoo multifasta_query

twofoo-clean:
	rm -fr twofoo twofoo_k31_r1 twofoo_k31_r1_hashval_k51/ \
        twofoo_k31_r1_multifasta/twofoo_k31_r1_search_oh0/
