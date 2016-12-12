all:

test:
	py.test

clean:
	rm -fr acido_bench/

sigs:
	sourmash compute -k 31 --dna data/acido-short.fa -o data/acido-short.fa.sig --name-from-first
	sourmash compute -k 31 --dna data/tr-1.fa -o data/tr-1.fa.sig --name-from-first
	sourmash compute -k 31 --dna data/tr-2.fa -o data/tr-2.fa.sig --name-from-first
