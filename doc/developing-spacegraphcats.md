# spacegraphcats: a brief developer guide

Spacegraphcats is written as a [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline.
The DAG below gives an overview of the steps in the pipeline and the scripts that execute them.

![](https://i.imgur.com/HQ21aM1.png)

## spacegraphcats files and what scripts produce them

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

