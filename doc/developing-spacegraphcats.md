# spacegraphcats: a brief developer guide

Spacegraphcats is written as a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline.
The DAG below gives an overview of the steps in the pipeline and the scripts that execute them.

![](https://i.imgur.com/HQ21aM1.png)

The spacegraphcats command line interface uses the [click API](https://click.palletsprojects.com/).
This also allows any snakemake flag to be added to the spacegraphcats CLI (e.g. `--unlock`).
 
## spacegraphcats files and what scripts produce them

* the `dory_k21/` directory contains the various cDBG files
    * `bcalm.inputlist.txt`
    * `bcalm.unitigs.fa`
    * `bcalm.unitigs.pickle`
    * `bcalm.unitigs.db` - database of unitigs in the cDBG. A FASTA file can be produced from this file by running `python -m spacegraphcats.cdbg.dump_contigs_db_to_fasta dory_k21/bcalm.unitigs.db`. 
    * `bcalm.unitigs.fa.sig` - sourmash signature for the cDBG nodes. Defaults to scaled = 1000. 
* the `dory_k21_r1/` directory contains the various files of the catlas constructed by spacegraphcats:
    * cdbg.gxt - the cDBG connection graph, in a custom format; produced by `bcalm_to_gxt.py`
    * reads.bgz.index - a BGZF file containing the cDBG unitigs; produced by `bcalm_to_gxt.py`
    * contigs.info.csv - summary information about the cDBG unitigs; produced by `bcalm_to_gxt.py`
    * contigs.indices - a numpy savez file containing mapping arrays; produced by `index_contigs_by_kmer.py`
    * contigs.mphf - Minimal Perfect Hash Function parameters for all of the k-mers in the cDBG; produced by `index_contigs_by_kmer.py`
    * contigs.sig - sourmash signature for cDBG nodes. Defaults to scaled = 1000. 
    * contigs.sizes - sizes of all cDBG nodes in pickle format; produced by `index_cdbg_by_kmer.py`
    * first_doms.txt - the dominating set information for the cDBG; produced by `spacegraphcats.catlas.catlas`
    * catlas.csv - the catlas for the cDBG; produced by `spacegraphcats.catlas.catlas`
    * commands.log - a partial log of all of the commands
* the `dory_k21_r1_search_oh0/` directory contains the output of a search:
    * results.csv - summary results for the queries (containment, similarity, etc.)
    * dory-head.fa.cdbg_ids.txt.gz - cDBG node IDs (unitig IDs) matching query
    * dory-head.fa.cdbg_ids.reads.gz - if `extract_reads` is used, reads from the original sequencing file that contain k-mers in the query neighborhood contigs 
    * dory-head.fa.contigs.sig - the sourmash signature of the entire match in the cDBG
    * dory-head.fa.frontier.txt.gz - catlas node IDs matching query
    * dory-head.fa.response.txt - response curve showing how much overhead is gained for each node
    * command.txt - a partial log of the commands run

## Spacegraphcats internal commands

For a brief look at some internal commands, see [this slideshow](https://hackmd.io/@ctb/BkDE7msPv).

### `extract_reads.py`: going from k-mer to cDBG/unitig ID to read

Spacegraphcats optionally outputs the reads that contain k-mers in a neighborhood. 
This requires being able to go from k-mer to cDBG unitig to read. 
Below we describe how is this implemented.
The code for that retrieval is in [extract_reads.py](https://github.com/spacegraphcats/spacegraphcats/blob/latest/spacegraphcats/search/extract_reads.py).

First, we put all of the reads in a BGZF file, which is a gzipped file that can be accessed by position. 
This lets us retrieve reads based solely on an offset into that file. 
See [make_bgzf.py](https://github.com/spacegraphcats/spacegraphcats/blob/latest/spacegraphcats/utils/make_bgzf.py) and [bgzf/](https://github.com/spacegraphcats/spacegraphcats/tree/latest/spacegraphcats/utils/bgzf) for our copy of the BGZF indexing code, taken from BioPython.

Then, we build a table that connects query k-mers to cDBG unitig IDs. 
Here the offset in the table is the MPHF of the query k-mer, built using bbhash, and the value in the table is the number of the cDBG unitig. 
This connects with internal spacegraphcats stuff.

Last but not least, we build a sqlite unitig-to-read offset table in [index_cdbg_by_kmer.py](https://github.com/spacegraphcats/spacegraphcats/blob/latest/spacegraphcats/cdbg/index_cdbg_by_kmer.py). 
This is a multimap table (read, unitig ID) that lets us query for the BGZF offset for all reads that belong to a given unitig ID. 
Once we have the read offsets for a given unitig ID, using the bgzf code we can reach into a sequence file and retrieve the relevant read(s) directly. 

## Test data sets

We use two test data sets, `dory` and `twofoo`. 
These data sets are packaged in the github repo under `tests/test-data`.
Dory is a subset of contigs from a Doryteuthis RNAseq assembly and is the smallest test set at 300KB. 
For a full example workflow using this data set, see the repository [spacegraphcats-dory-example](https://github.com/spacegraphcats/spacegraphcats-dory-example).

Twofoo is a subset of reads from a synthetic metagenome sequenced by [Shakya et al., 2013](https://www.ncbi.nlm.nih.gov/pubmed/23387867).
The metagenome was constructed by aliquoting isolate DNA into a mixture and sequencing this mixture.
We often refer to this metagenom as "podar", the last name of the senior author on the Shakya et al. 2013 paper.
The test data is composed of two sets of reads, `shew-reads` and `akker-reads`.
`shew-reads.abundtrim.gz` is a collection of reads from podar data that maps to the Shewanella OS223 genome via `bwa aln`.  
Note that there is significant overlap with the *Shewanella baltica* OS185 genome; meaning these reads have significant strain variation.
`akker-reads.abundtrim.gz` is a collection of reads from podar that maps to *Akkermansia muciniphila* ATCC BAA-835 via `bwa aln`.
The `akker-reads` do not have strain variation.
For a full example workflow using this data set, see the repository [spacegraphcats-twofoo-example](https://github.com/spacegraphcats/spacegraphcats-twofoo-example).
