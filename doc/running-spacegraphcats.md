# spacegraphcats: a run guide

## Installing the spacegraphcats software and its dependencies

### If starting on a blank Ubuntu machine

e.g. on AWS, ubuntu/images/hvm-ssd/ubuntu-xenial-16.04-amd64-server-20180126 (ami-79873901), you'll need to make sure you have Python 3, a dev environment, and other stuff:

```
sudo apt-get update
sudo apt-get -y install python3 python3-dev zlib1g-dev g++ \
    python3-venv make cmake
```

and then do

```
python3 -m venv catsenv
```

instead of the first command below.

### First, clone repo and configure/install requirements

Change to a working directory, and create a virtualenv; you'll need Python 3.5 or up.

```
python -m virtualenv -p python3.5 catsenv 
```

Activate the virtualenv:
```
. catsenv/bin/activate
```

Next, clone spacegraphcats.
```
git clone https://github.com/spacegraphcats/spacegraphcats/
```

and now install the requirements:

```
cd spacegraphcats
pip install -U setuptools pip
pip install Cython
pip install -r requirements.txt
```

You will also need to install BCALM and put the bcalm binary in your path; [see instructions](https://github.com/GATB/bcalm#installation).

### Now, run a small test.

```
conf/run dory-test search
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

* the `dory/` directory contains the BCALM assembly of the input files into a compact De Bruijn graph; the key file here is `dory/bcalm.dory.k21.unitigs.fa`. There is also an output log file, `bcalm.dory.k21.unitigs.fa.log.txt`, that contains the console output of BCALM's run.
* the `dory_k21_r1/` directory contains the various files of the catlas constructed by spacegraphcats:
    * cdbg.gxt
    * contigs.fa.gz
    * contigs.fa.gz.indices
    * contigs.fa.gz.info.csv
    * contigs.fa.gz.mphf
    * first_doms.txt
    * catlas.csv
    * commands.log
* the `dory_k21_r1_search_oh0/` directory contains the output of a search:
    * results.csv
    * dory-head.fa.cdbg_ids.txt.gz
    * dory-head.fa.contigs.sig
    * dory-head.fa.frontier.txt.gz
    * dory-head.fa.response.txt
    * command.txt

## Configuring and running spacegraphcats itself

There are three important top-level files.

1. The script `conf/run` uses [snakemake](https://snakemake.readthedocs.io/en/stable/) to run the pipeline.

2. The Snakefile itself is under `conf/Snakefile`, and is not intended to be particularly human readable :). The snakemake run is configured by a JSON config file.

3. The configuration files are relatively simple, and are (by convention) contained in the `conf/` directory, although you can put them elsewhere.

To run spacegraphcats on a new data set, you will need to write a new config file.

## Config files: the `dory` example

Config files look like this:
```
{
    "catlas_base": "dory",
    "input_sequences": ["data/dory-subset.fa"],
    "ksize": 21,
    "radius": 1,

    "searchquick": ['data/dory-head.fa'],
    "searchseeds": "43",
    "search": ['data/dory-head.fa'],
}
```

Here, the `catlas_base`, `ksize`, and `radius` are used to name output directories - hence `dory/` and `dory_k21_r1/`, above. `catlas_base` is arbitrary, `ksize` is the k-mer size (we suggest 21 or 31), and `radius` is the domset radius (we suggest 1).

The `input_sequences` contains the set of input files (FASTA or FASTQ, potentially gzipped) to use to build the compact De Bruijn graph. Note, these should be error-trimmed, or otherwise you will have very big cDBGs...

The `searchquick` and `search` specify the *query files* that you are using to search the catlas.

So, in brief, this config file:

* builds a catlas with k=21 and r=1 from `data/dory-subset.fa`;
* searches with `data/dory-head.fa` when you call either `search` or `searchquick`.

## A larger example: `twofoo`

Here you can see another config file, for a synthetic mixture
of two genomes (Akkermansia and Shewanella baltica OS 223) from the Shakya et al., 2014, benchmark study:
```
{
    "catlas_base": "twofoo",
    "input_sequences": ["twofoo.fq.gz"],
    "ksize": 31,
    "radius": 1,

    "searchquick": ["data/63.fa.gz"],

    "search": ["data/2.fa.gz", "data/47.fa.gz", "data/63.fa.gz"],

}
```

To run this, you'll first need to download and prepare the input file, `twofoo.fq.gz`:

```
make twofoo.fq.gz
```

which will take a few minutes.

Then, run:
```
conf/run twofoo search
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

You can get the cDBG sequences by using `search.extract_contigs`; this will extract all of the contigs matched by the `2.fa.gz` query from the `twofoo_k31_r1` catlas:

```
python -m search.extract_contigs twofoo_k31_r1 twofoo_k31_r1_search_oh0/2.fa.gz.cdbg_ids.txt.gz -o xxx.2.query.fa
```

You can also extract the reads, but this is a little rougher at the moment. For that, you will first need to put all the reads into a .bgz file:

```
python -m search.make_bgzf twofoo.fq.gz
```

and then index that file by contig:

```
python -m search.label_cdbg twofoo_k31_r1 twofoo.fq.gz.bgz twofoo.labels
```

and NOW FINALLY you can run
```
python -m search.extract_reads twofoo.fq.gz.bgz twofoo.labels twofoo_k31_r1_search_oh0/2.fa.gz.cdbg_ids.txt.gz
```
to extract the reads belonging to the cDBG ids returned from the query with `data/2.fa.gz`.

## Other information

### The `run` script

The `run` script has several targets, in addition to `search` and `searchquick`.

* `conf/run twofoo build` will build the catlas
* `conf/run twofoo clean` should remove the build targets.

You can also specify `--radius <n>` to override the radius defined in the JSON config file; `--overhead <fraction>` to specify an overhead for searches; and `--experiment foo` to append a `_foo` to the search directory.)

Last, but not least: snakemake locks the directory to make sure processes don't step on each other. This is important when catlases need to be built (you don't want two different `search` commands stepping on each other during the catlas building phase) but once you have built catlases you can do searches in parallel.  To enable this add `--nolock` to the run command. 

## Characterizing the catlas

### Extracting high articulated bits of the cDBG

First, build the twofoo data set with r5:

```
conf/run twofoo build --radius=5
```

Then, extract nodes with many cDBG nodes and few k-mers (by ratio):

```
python -m search.extract_nodes_by_shadow_ratio twofoo_k31_r5 zzz.fq
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
