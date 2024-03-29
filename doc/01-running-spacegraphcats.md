# spacegraphcats: a run guide

## Installing the spacegraphcats software and its dependencies

Please see [Installing spacegraphcats](00-installing-spacegraphcats.md).

## Running spacegraphcats search & output files

```
python -m spacegraphcats run dory-test search
```

This should run in a few seconds, and you should see something like this in the output:

```
...
query inclusion by retrieved contigs: 100.000%
query similarity to retrieved contigs: 100.000%
total time: 0.0s
Finished job 0.
6 of 6 steps (100%) done
```

You will have a bunch of new output files:

* the `dory_k21/` directory contains the BCALM assembly of the input files into a compact De Bruijn graph; the key file here is `dory_k21/bcalm.unitigs.db`. There is also an output log file, `bcalm.log.txt`, that contains the console output of BCALM's run.
* the `dory_k21_r1/` directory contains the various files of the catlas constructed by spacegraphcats:
    * cdbg.gxt - the cDBG connection graph, in a custom format; produced by `bcalm_to_gxt.py`
    * contigs.info.csv - summary information about the cDBG unitigs; produced by `bcalm_to_gxt.py`
    * contigs.indices - a numpy savez file containing mapping arrays; produced by `index_contigs_by_kmer.py`
    * contigs.mphf - Minimal Perfect Hash Function parameters for all of the k-mers in the cDBG; produced by `index_contigs_by_kmer.py`
    * contigs.sig - sourmash signature for cDBG nodes. Defaults to scaled = 1000.
    * contigs.sizes - sizes of all cDBG nodes in pickle format; produced by `index_cdbg_by_kmer.py`
    * first_doms.txt - the dominating set information for the cDBG; produced by `spacegraphcats/catlas/catlas.py`
    * catlas.csv - the catlas for the cDBG; produced by `spacegraphcats/catlas/catlas.py`
    * commands.log - a partial log of all of the commands
* the `dory_k21_r1_search_oh0/` directory contains the output of a search:
    * results.csv - summary results for the queries (containment, similarity, etc.)
    * dory-head.fa.cdbg_ids.txt.gz - cDBG node IDs (unitig IDs) matching query
    * dory-head.fa.cdbg_ids.reads.gz - if `extract_reads` is used, reads from the original sequencing file that contain k-mers in the query neighborhood contigs. Reads will be output in the same format as the input data, e.g. in FASTQ if the original data was FASTQ, or FASTA if the original data was FASTA.
    * dory-head.fa.contigs.sig - the sourmash signature of the entire match in the cDBG
    * dory-head.fa.frontier.txt.gz - catlas node IDs matching query
    * dory-head.fa.response.txt - response curve showing how much overhead is gained for each node
    * command.txt - a partial log of the commands run

The `results.csv` file contains the following information:

+ `query`:  the file path to the query used by spacegraphcats 
+ `containment`: the containment of the query in the retrieved neighborhood
+ `similarity`: the similarity of the query to the retrieved neighborhood
+ `bp`: the number of basepairs in the neighborhood
+ `contigs`: the number of unitigs contained in the query from the cDBG
+ `ksize`: the k-mer size used to perform the query
+ `num_query_kmers`: the number of k-mers in the query 
+ `best_containment`: tbd
+ `cdbg_min_overhead`: tbd
+ `catlas_min_overhead`: tbd

Note that:

+ if the query Q is entirely in the indexed metagenome M, then 1/similarity is the relative neighborhood size, because |Q intersect M| is |Q|, and |Q union M| is |M|, so 1/similarity is |M|/|Q|. This does not hold if the query is not entirely within the indexed metagenome.
+ `num_query_kmers` is |Q|, not |Q intersect M|

## Configuring and running spacegraphcats itself

There are three important top-level files.

1. The script `python -m spacegraphcats` uses [snakemake](https://snakemake.readthedocs.io/en/stable/) to run the pipeline.

2. The Snakefile itself is under `conf/Snakefile`, and is not intended to be particularly human readable :). The snakemake run is configured by a [YAML config file](https://www.cloudbees.com/blog/yaml-tutorial-everything-you-need-get-started/).

3. The configuration files are relatively simple, and are (by convention) contained in the `conf/` directory, although you can put them elsewhere. The filenames should end in either `.yaml` or `.conf`.

To run spacegraphcats on a new data set, you will need to write a new
config file. An empty template is available as
`spacegraphcats/conf/empty.yaml`, and you can display the empty config
file with `python -m spacegraphcats showconf empty`.

## Config files: the `dory` example

Config files look like this:
```
catlas_base: dory
input_sequences:
- data/dory-subset.fq
keep_graph_pendants: true
paired_reads: false
ksize: 21
radius: 1
search:
- data/dory-head.fa
```

Here, the `catlas_base`, `ksize`, and `radius` are used to name output directories - hence `dory/` and `dory_k21_r1/`, above. `catlas_base` is arbitrary, `ksize` is the k-mer size (we suggest 21 or 31), and `radius` is the domset radius (we suggest 1).

The `input_sequences` contains the set of input files (FASTA or FASTQ, potentially gzipped) to use to build the compact De Bruijn graph. Note, these should be k-mer-error-trimmed, or otherwise you will have very big cDBGs... see [the error trimming information from Brown et al. for an example](https://github.com/spacegraphcats/2018-paper-spacegraphcats/tree/master/pipeline-base#pipeline-to-generate-files-needed-for-figures-2-and-3). Note that the `trim-low-abund.py` command is in the khmer package, which is a requirement for spacegraphcats, so you will already have it installed.

The `search` entry specifies the *query files* that you are using to search the catlas.

So, in brief, this config file:

* builds a catlas with k=21 and r=1 from `data/dory-subset.fa`;
* searches with `data/dory-head.fa` when you call `search`.

## A larger example: `twofoo`

Here you can see another config file, for a synthetic mixture
of two genomes (Akkermansia and Shewanella baltica OS 223) from the Shakya et al., 2014, benchmark study:
```
catlas_base: twofoo
input_sequences:
- twofoo.fq.gz
ksize: 31
radius: 1
search:
- data/2.fa.gz
- data/47.fa.gz
- data/63.fa.gz
```

To run this, you'll first need to download and prepare the input file, `twofoo.fq.gz`:

```
make twofoo.fq.gz
```

which will take a few minutes.

Then, run:
```
python -m spacegraphcats run twofoo search
```

which will generate searches of the twofoo synthetic data set with `data/2.fa.gz`, `data/47.fa.gz`, and `data/63.fa.gz`.

It should take about 5 minutes on a relatively recent laptop, and will require ~1 GB of RAM.

## Interpreting the results of a run

Start up a Jupyter notebook, and run this in a cell:

```
%matplotlib inline
import pandas
from matplotlib import pyplot

results = pandas.read_csv('twofoo_k31_r1_search_oh0/results.csv')

fig, axes = pyplot.subplots(nrows=1, ncols=2, figsize=(6, 3))
axes[0].violinplot(results.similarity)
axes[1].violinplot(results.containment)

axes[0].axis(ymin=0.0)
axes[0].axis(ymax=1.0)

axes[0].get_xaxis().set_ticks([])
axes[1].get_xaxis().set_ticks([])
axes[1].get_yaxis().set_ticks([])
    

axes[0].set_xlabel('similarity')
axes[1].set_xlabel('containment')
```
This will plot the distribution of similarity and containment in the results.

## Extracting the sequences that match search results

You can get the cDBG unitigs for the searches, and the reads
corresponding to them, by using the targets `extract_reads` and
`extract_contigs`.

```
python -m spacegraphcats run twofoo extract_contigs extract_reads
```

This will produce the files:

```
twofoo_k31_r1_search_oh0/2.fa.gz.cdbg_ids.contigs.fa.gz
twofoo_k31_r1_search_oh0/47.fa.gz.cdbg_ids.contigs.fa.gz
twofoo_k31_r1_search_oh0/63.fa.gz.cdbg_ids.contigs.fa.gz

twofoo_k31_r1_search_oh0/2.fa.gz.cdbg_ids.reads.gz
twofoo_k31_r1_search_oh0/47.fa.gz.cdbg_ids.reads.gz
twofoo_k31_r1_search_oh0/63.fa.gz.cdbg_ids.reads.gz
```
which are (respectively) the contigs for the neighborhoods around each
query, and the reads for the neighborhoods around each query.

## Other information

### The `spacegraphcats` script

The `spacegraphcats` script has several targets, in addition to `search`.

* `python -m spacegraphcats run twofoo build` will build the catlas
* `python -m spacegraphcats run twofoo clean` should remove the build targets.
* `python -m spacegraphcats run twofoo extract_contigs` -- get contigs for search results; see above.
* `python -m spacegraphcats run twofoo extract_reads` -- get reads for search results; see above.

You can also specify `--radius <n>` to override the radius defined in the YAML config file and `--experiment foo` to append a `_foo` to the search directory.)

Last, but not least: snakemake locks the directory to make sure processes don't step on each other. This is important when catlases need to be built (you don't want two different `search` commands stepping on each other during the catlas building phase) but once you have built catlases you can do searches in parallel.  To enable this add `--nolock` to the run command. 

## Characterizing the catlas

(The below may not be working. - CTB 7/22/2018.)

### Extracting high articulated bits of the cDBG

First, build the twofoo data set with r5:

```
python -m spacegraphcats run twofoo build --radius=5
```

Then, extract nodes with many cDBG nodes and few k-mers (by ratio):

```
python -m spacegraphcats.search.extract_nodes_by_shadow_ratio twofoo_k31_r5 zzz.fq
```

Now, look at the content of the extracted nodes -- the presence of the Akkermansia genome is essentially nil,
because there is no strain variation in this part of the graph / the assembly is quite good:

```
sourmash search --containment zzz.fq.sig data/2-akker.sig --threshold=0
```

should yield
```
similarity   match
----------   -----
  0.2%       CP001071.1 Akkermansia muciniphila ATCC BAA-835, complete...
```

while the two Shewanella genomes interfere with each other, creating a highly articulated cDBG that leads to poor assembly.  Thus the command above extracts the Shewanella bits of the catlas preferentially; so

```
sourmash search --containment zzz.fq.sig data/47-os185.sig --threshold=0
```

should yield

```
similarity   match
----------   -----
 54.1%       NC_009665.1 Shewanella baltica OS185, complete genome
```

and

```
sourmash search --containment zzz.fq.sig data/63-os223.sig --threshold=0
```

should yield

```
similarity   match
----------   -----
 59.1%       NC_011663.1 Shewanella baltica OS223, complete genome
```
