all:
	python setup.py build_ext -i

test:
	python setup.py build_ext -i
	py.test

clean:
	rm -fr acido_bench/

bench: acido_bench_target

15genome_target: 15genome/results_bestnode.csv \
	15genome/results_searchlevel.5.csv \
	15genome/results_gathermins.5.csv \
	15genome/results_gathermins2.5.csv

15genome/15genome.gxt: data/15genome.fa.gz
	./walk-dbg.py -x 8e9 data/15genome.fa.gz --label-linear-segments --label

15genome/15genome.catlas.5.gxt: 15genome/15genome.gxt
	./build-catlas.py 15genome 5

15genome/results_bestnode.csv: 15genome/15genome.catlas.5.gxt
	rm -f 15genome/results_bestnode.csv
	search-for-domgraph-nodes-multi.py --quiet 15genome 5 ./data/15genome.fa.{?,??}.sig.dump.txt --append-csv 15genome/results_bestnode.csv

15genome/results_searchlevel.5.csv: 15genome/15genome.catlas.5.gxt
	rm -f 15genome/results_searchlevel*.csv
	search-for-domgraph-nodes-multi.py --quiet 15genome 5 ./data/15genome.fa.{?,??}.sig.dump.txt --append-csv 15genome/results_searchlevel.3.csv --strategy searchlevel --searchlevel 3
	search-for-domgraph-nodes-multi.py --quiet 15genome 5 ./data/15genome.fa.{?,??}.sig.dump.txt --append-csv 15genome/results_searchlevel.5.csv --strategy searchlevel --searchlevel 5

15genome/results_gathermins.5.csv: 15genome/15genome.catlas.5.gxt
	rm -f 15genome/results_gathermins*.csv
	search-for-domgraph-nodes-multi.py --quiet 15genome 5 ./data/15genome.fa.{?,??}.sig.dump.txt --append-csv 15genome/results_gathermins.3.csv --strategy gathermins --searchlevel 3
	search-for-domgraph-nodes-multi.py --quiet 15genome 5 ./data/15genome.fa.{?,??}.sig.dump.txt --append-csv 15genome/results_gathermins.5.csv --strategy gathermins --searchlevel 5

15genome/results_gathermins2.5.csv: 15genome/15genome.catlas.5.gxt
	rm -f 15genome/results_gathermins2*.csv
	search-for-domgraph-nodes-multi.py --quiet 15genome 5 ./data/15genome.fa.{?,??}.sig.dump.txt --append-csv 15genome/results_gathermins2.3.csv --strategy gathermins2 --searchlevel 3
	search-for-domgraph-nodes-multi.py --quiet 15genome 5 ./data/15genome.fa.{?,??}.sig.dump.txt --append-csv 15genome/results_gathermins2.5.csv --strategy gathermins2 --searchlevel 5

acido_bench_target: acido_bench/results_bestnode.csv \
	acido_bench/results_searchlevel.csv \
	acido_bench/results_gathermins.csv \
	acido_bench/results_gathermins2.csv

acido_bench/acido_bench.gxt:
	-rm -fr acido_bench
	./walk-dbg.py --label -x 1e9 -k 31 data/acido-chunk*.fa.gz -o acido_bench --label-linear-segments

acido_bench/acido_bench.catlas.5.gxt: acido_bench/acido_bench.gxt
	./build-catlas.py acido_bench 5

acido_bench/results_bestnode.csv: acido_bench/acido_bench.catlas.5.gxt
	rm -f acido_bench/results_bestnode.csv
	search-for-domgraph-nodes-multi.py --quiet acido_bench 5 ./data/acido-chunk*.fa.sig.dump.txt --append-csv acido_bench/results_bestnode.csv

acido_bench/results_searchlevel.csv: acido_bench/acido_bench.catlas.5.gxt
	rm -f acido_bench/results_searchlevel.csv
	search-for-domgraph-nodes-multi.py --quiet acido_bench 5 ./data/acido-chunk*.fa.sig.dump.txt --append-csv acido_bench/results_searchlevel.csv --strategy searchlevel --searchlevel 3

acido_bench/results_gathermins.csv: acido_bench/acido_bench.catlas.5.gxt
	rm -f acido_bench/results_gathermins.csv
	search-for-domgraph-nodes-multi.py --quiet acido_bench 5 ./data/acido-chunk*.fa.sig.dump.txt --append-csv acido_bench/results_gathermins.csv --strategy gathermins --searchlevel 3

acido_bench/results_gathermins2.csv: acido_bench/acido_bench.catlas.5.gxt
	rm -f acido_bench/results_gathermins2.csv
	search-for-domgraph-nodes-multi.py --quiet acido_bench 5 ./data/acido-chunk*.fa.sig.dump.txt --append-csv acido_bench/results_gathermins2.csv --strategy gathermins2 --searchlevel 3
