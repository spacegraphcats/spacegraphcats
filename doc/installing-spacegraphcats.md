# spacegraphcats: installation

## First, install the dependencies.

spacegraphcats relies on several dependencies - Python, a C/C++ environment,
and bcalm 2, in particular.

There are several ways to install these dependencies. We recommend
using conda, but include instructions for using a Python virtual
environment or a blank Ubuntu machine as well.

### 1. Installing dependencies using conda

If you have [conda](https://docs.conda.io/en/latest/) installed, you can
use packages from conda-forge and [bioconda](https://docs.conda.io/en/latest/)
to install all the dependencies for spacegraphcats.

First, create a new 'sgc' environment, activate, and install the necessary
packages:

```
conda create -n sgc python==3@@
conda activate sgc
conda install -c conda-forge -c bioconda ... bcalm
``

Then, follow the instructions below for installing spacegraphcats.

## 2. Installing dependencies in a virtual environment

If you already have a functioning Python >= 3.5 along with a C/C++
development environment, you can install the dependencies in a
venv.

Change to a working directory, and create a virtualenv:

```
python -m virtualenv -p python3.5 catsenv
```

Activate the virtualenv, upgrade pip, and install Cython:
```
. catsenv/bin/activate
pip install -U setuptools pip
pip install Cython
```

**Note:** You will also need to install bcalm 2; please following
[their install instructions](https://github.com/GATB/bcalm#installation).

Now, follow the instructions below for installing spacegraphcats.

## 3. Installing dependencies on a blank Ubuntu machine

If you're starting e.g. on a blank AWS instance such as
ubuntu/images/hvm-ssd/ubuntu-xenial-16.04-amd64-server-20180126
(ami-79873901), you'll need to make sure you have Python 3, a dev
environment, and other stuff:

```
sudo apt-get update
sudo apt-get -y install python3 python3-dev zlib1g-dev g++ \
    python3-venv make cmake
```

Now create a virtualenv named catsenv and activate it; then upgrade
and install a few things.

```
python3 -m venv catsenv
. catsenv/bin/activate
pip install -U setuptools pip
pip install Cython
```

**Note:** You will also need to install bcalm 2; please following
[their install instructions](https://github.com/GATB/bcalm#installation).

Finally, follow the instructions below for installing spacegraphcats.

## After installing dependencies, install spacegraphcats

Once you have the basic dependencies, go ahead and clone spacegraphcats:
```
git clone https://github.com/spacegraphcats/spacegraphcats/
```

and install the requirements:

```
cd spacegraphcats
pip install -r requirements.txt
```

This will take a few minutes.

## Test spacegraphcats by running a small test: dory.

```
python -m spacegraphcats dory-test search
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
    
## Run a bigger test: twofoo.

Now let's check the full pipeline for a synthetic mixture of two
genomes (Akkermansia and Shewanella baltica OS 223) from the Shakya et
al., 2014, benchmark study.

To run this, you'll first need to download and prepare the input file, `twofoo.fq.gz`:

```
make twofoo.fq.gz
```

which will take a few minutes.

Then, run:
```
python -m spacegraphcats twofoo search
```

which will generate searches of the twofoo synthetic data set with `data/2.fa.gz`, `data/47.fa.gz`, and `data/63.fa.gz`.

It should take about 5 minutes on a relatively recent laptop, and will require ~1 GB of RAM.

If this succeeds, your install is all good - give yourself a high five!
you've made it!

See [running spacegraphcats](running-spacegraphcats.md) for next steps and
more information on config files.
