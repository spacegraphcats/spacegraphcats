# spacegraphcats: installation

## If starting on a blank Ubuntu machine

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

## First, clone repo and configure/install requirements

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
pip install https://github.com/dib-lab/khmer/archive/master.zip
```

You will also need to install BCALM and put the bcalm binary in your path; [see instructions](https://github.com/GATB/bcalm#installation).

### Now, run a small test: dory.

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
conf/run twofoo search
```

which will generate searches of the twofoo synthetic data set with `data/2.fa.gz`, `data/47.fa.gz`, and `data/63.fa.gz`.

It should take about 5 minutes on a relatively recent laptop, and will require ~1 GB of RAM.

If this succeeds, your install is all good - give yourself a high five!
you've made it!

See [running spacegraphcats](running-spacegraphcats.md) for next steps and
more information on config files.
