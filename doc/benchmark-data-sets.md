# Benchmark data sets

## 15genome

The '15genome' data set contains the first 15 genomes used to build
the mock community in the
[Shakya et al., 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3665634/)
paper.  This is a labeled data set - the cDBG components corresponding to
each genome are labeled, and MinHash signatures are provided for use.
k
CTB TODO: provide full build instructions, including signature computation.

The catlas, computed with r=3, is available for download
[here](http://spacegraphcats.ucdavis.edu.s3.amazonaws.com/15genome.3.catlas.tar.gz).  To download and unpack it:

    curl -O http://spacegraphcats.ucdavis.edu.s3.amazonaws.com/15genome.3.catlas.tar.gz
    tar xzf 15genome.3.catlas.tar.gz
    
This creates the `15genome` directory.  The signatures are available within
the directory as `15genome/15genome.fa.gz.sig`.

To run against the entire set of signatures, do:

    ./search-for-domgraph-nodes.py 15genome 3 15genome/15genome.fa.gz.sig
    
Add `--append-csv out.csv` to get the output in spreadsheet form.

## Appendices

### Appendix: Installing the requirements

Starting from a blank Ubuntu 15.10 install, add some packages:

    sudo apt-get update && sudo apt-get -y install python3.5-dev \
       python3-virtualenv python3-matplotlib python3-numpy g++ make
       
Create a new virtualenv:

    cd
    python3.5 -m virtualenv env -p python3.5 --system-site-packages
    . env/bin/activate

Install many things:

    pip install screed pytest PyYAML
    pip install git+https://github.com/dib-lab/khmer.git
    pip install git+https://github.com/dib-lab/sourmash.git

Now clone and install spacegraphcats:

    git clone https://github.com/spacegraphcats/spacegraphcats

and run the tests to be sure:

    cd spacegraphcats
    make test
    
You should be done!

## Appendix: Updating the software

The two most likely packages for updates will be khmer and sourmash.  To
update, run these commands from within the virtualenv:

    pip install git+https://github.com/dib-lab/khmer.git
    pip install git+https://github.com/dib-lab/sourmash.git
