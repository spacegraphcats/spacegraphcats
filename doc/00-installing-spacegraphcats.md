# spacegraphcats: installation

## Step I: install the dependencies.

spacegraphcats relies on several dependencies - Python, a C/C++ environment,
and bcalm 2, in particular.

There are several ways to install these dependencies. We recommend
using conda, but include instructions for using a Python virtual
environment or a blank Ubuntu machine as well.

### 1. Installing dependencies using conda

If you have [conda](https://docs.conda.io/en/latest/) installed, you can
use packages from conda-forge and [bioconda](https://docs.conda.io/en/latest/)
to install all the dependencies for spacegraphcats.

Start by cloning the spacegraphcats repository:

```
git clone https://github.com/spacegraphcats/spacegraphcats/
```

Now, create a new 'sgc' environment with the necessary packages
packages:

```
conda env create -f spacegraphcats/environment.yml -n sgc
```

Once it's created, activate the environment you just created:

```
conda activate sgc
```

and install spacegraphcats from the development directory:
```
pip install -e ./spacegraphcats/
```

and voila, done!

## 2. Installing dependencies in a virtual environment

If you already have a functioning Python >= 3.7 along with a C/C++
development environment, you can install the dependencies in a
venv.

Change to a working directory, and create a virtualenv:

```
python -m virtualenv -p python3.7 catsenv
```

Activate the virtualenv, upgrade pip, and install Cython:
```
. catsenv/bin/activate
pip install -U setuptools pip
pip install Cython
```

**Note:** You will also need to install bcalm 2; please following
[their install instructions](https://github.com/GATB/bcalm#installation).

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

Finally, install the git repo in developer mode:

```
pip install -e ./spacegraphcats/
```

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


Finally, install the git repo in developer mode:

```
pip install -e ./spacegraphcats/
```

## Step II: Test spacegraphcats by running a small test: dory.

In the `spacegraphcats/` top level directory (containing e.g. `README.md`),
run:

```
python -m spacegraphcats run dory-test search
```

This should run in a few seconds, and you should see something like this in the output:

```
...
=> containment: 0.0%
[Fri Aug 14 07:53:04 2020]
Finished job 0.
4 of 4 steps (100%) done
Complete log: /Users/t/dev/spacegraphcats/.snakemake/log/2020-08-14T075302.093387.snakemake.log

-------- DONE --------

catlas output directory: dory_k21_r1
search output directory: dory_k21_r1_search_oh0
```

You will have a bunch of new output files:

* the `dory_k21/` directory contains the BCALM assembly of the input files into a compact De Bruijn graph; the key file here is `dory_k21/bcalm.unitigs.db`. There is also an output log file, `bcalm.log.txt`, that contains the console output of BCALM's run.
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
    
## Step III: Run a bigger test: twofoo.

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
python -m spacegraphcats run twofoo search
```

which will generate searches of the twofoo synthetic data set with `data/2.fa.gz`, `data/47.fa.gz`, and `data/63.fa.gz`.

It should take about 5 minutes on a relatively recent laptop, and will require ~1 GB of RAM.

If this succeeds, your install is all good - give yourself a high five!
you've made it!

See [running spacegraphcats](01-running-spacegraphcats.md) for next steps and
more information on config files.
