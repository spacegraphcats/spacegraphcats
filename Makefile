all:

test:
	py.test

clean:
	rm -fr acido_bench/

bench: acido_bench_target

acido_bench_target: acido_bench/results_bestnode.csv \
	acido_bench/results_searchlevel.csv \
	acido_bench/results_gathermins.csv \
	acido_bench/results_gathermins2.csv

acido_bench/acido_bench.gxt:
	-rm -fr acido_bench
	./walk-dbg.py --label -x 1e9 -k 31 data/acido-chunk*.fa.gz -o acido_bench

acido_bench/acido_bench.catlas.5.gxt: acido_bench/acido_bench.gxt
	./build-catlas.py acido_bench 5

acido_bench/results_bestnode.csv: acido_bench/acido_bench.catlas.5.gxt
	rm -f acido_bench/results_bestnode.csv
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk1.fa.sig.dump.txt 1 --append-csv acido_bench/results_bestnode.csv
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk2.fa.sig.dump.txt 2 --append-csv acido_bench/results_bestnode.csv
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk3.fa.sig.dump.txt 3 --append-csv acido_bench/results_bestnode.csv
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk4.fa.sig.dump.txt 4 --append-csv acido_bench/results_bestnode.csv
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk5.fa.sig.dump.txt 5 --append-csv acido_bench/results_bestnode.csv
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk6.fa.sig.dump.txt 6 --append-csv acido_bench/results_bestnode.csv
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk7.fa.sig.dump.txt 7 --append-csv acido_bench/results_bestnode.csv
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk8.fa.sig.dump.txt 8 --append-csv acido_bench/results_bestnode.csv

acido_bench/results_searchlevel.csv: acido_bench/acido_bench.catlas.5.gxt
	rm -f acido_bench/results_searchlevel.csv
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk1.fa.sig.dump.txt 1 --append-csv acido_bench/results_searchlevel.csv --strategy searchlevel --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk2.fa.sig.dump.txt 2 --append-csv acido_bench/results_searchlevel.csv --strategy searchlevel --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk3.fa.sig.dump.txt 3 --append-csv acido_bench/results_searchlevel.csv --strategy searchlevel --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk4.fa.sig.dump.txt 4 --append-csv acido_bench/results_searchlevel.csv --strategy searchlevel --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk5.fa.sig.dump.txt 5 --append-csv acido_bench/results_searchlevel.csv --strategy searchlevel --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk6.fa.sig.dump.txt 6 --append-csv acido_bench/results_searchlevel.csv --strategy searchlevel --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk7.fa.sig.dump.txt 7 --append-csv acido_bench/results_searchlevel.csv --strategy searchlevel --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk8.fa.sig.dump.txt 8 --append-csv acido_bench/results_searchlevel.csv --strategy searchlevel --searchlevel 3

acido_bench/results_gathermins.csv: acido_bench/acido_bench.catlas.5.gxt
	rm -f acido_bench/results_gathermins.csv
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk1.fa.sig.dump.txt 1 --append-csv acido_bench/results_gathermins.csv --strategy gathermins --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk2.fa.sig.dump.txt 2 --append-csv acido_bench/results_gathermins.csv --strategy gathermins --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk3.fa.sig.dump.txt 3 --append-csv acido_bench/results_gathermins.csv --strategy gathermins --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk4.fa.sig.dump.txt 4 --append-csv acido_bench/results_gathermins.csv --strategy gathermins --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk5.fa.sig.dump.txt 5 --append-csv acido_bench/results_gathermins.csv --strategy gathermins --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk6.fa.sig.dump.txt 6 --append-csv acido_bench/results_gathermins.csv --strategy gathermins --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk7.fa.sig.dump.txt 7 --append-csv acido_bench/results_gathermins.csv --strategy gathermins --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk8.fa.sig.dump.txt 8 --append-csv acido_bench/results_gathermins.csv --strategy gathermins --searchlevel 3

acido_bench/results_gathermins2.csv: acido_bench/acido_bench.catlas.5.gxt
	rm -f acido_bench/results_gathermins2.csv
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk1.fa.sig.dump.txt 1 --append-csv acido_bench/results_gathermins2.csv --strategy gathermins2 --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk2.fa.sig.dump.txt 2 --append-csv acido_bench/results_gathermins2.csv --strategy gathermins2 --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk3.fa.sig.dump.txt 3 --append-csv acido_bench/results_gathermins2.csv --strategy gathermins2 --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk4.fa.sig.dump.txt 4 --append-csv acido_bench/results_gathermins2.csv --strategy gathermins2 --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk5.fa.sig.dump.txt 5 --append-csv acido_bench/results_gathermins2.csv --strategy gathermins2 --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk6.fa.sig.dump.txt 6 --append-csv acido_bench/results_gathermins2.csv --strategy gathermins2 --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk7.fa.sig.dump.txt 7 --append-csv acido_bench/results_gathermins2.csv --strategy gathermins2 --searchlevel 3
	search-for-domgraph-nodes.py --quiet acido_bench 5 ./data/acido-chunk8.fa.sig.dump.txt 8 --append-csv acido_bench/results_gathermins2.csv --strategy gathermins2 --searchlevel 3
