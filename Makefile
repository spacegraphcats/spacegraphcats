all: test

flakes:
	pyflakes search/*.py

lint:
	mypy spacegraphcats/*.py search/*.py --ignore-missing-imports

test:
	py.test spacegraphcats search

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
	conf/run dory-test

twofoo-test:
	python -m search.extract_reads_by_shadow_ratio twofoo twofoo.fq.gz.bgz twofoo.labels twofoo.shadow.out.fa -k 31
	python -m search.extract_nodes_by_query twofoo foo --query data/2.fa.gz -k 31
