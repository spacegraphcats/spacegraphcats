# spacegraphcats: a brief developer guide

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
