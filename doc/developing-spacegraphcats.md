# spacegraphcats: a brief developer guide

## spacegraphcats files and what scripts produce them

* the `dory_k21_r1/` directory contains the various files of the catlas constructed by spacegraphcats:
    * cdbg.gxt - the cDBG connection graph, in a custom format; produced by `bcalm_to_gxt.py`
    * contigs.fa.gz - a BGZF file containing the cDBG unitigs; produced by `bcalm_to_gxt.py`
    * contigs.fa.gz.info.csv - summary information about the cDBG unitigs; produced by `bcalm_to_gxt.py`
    * contigs.fa.gz.indices - a numpy savez file containing mapping arrays; produced by `index_contigs_by_kmer.py`
    * contigs.fa.gz.mphf - Minimum Perfect Hash Function parameters for all of the k-mers in the cDBG; produced by `index_contigs_by_kmer.py`
    * first_doms.txt - the dominating set information for the cDBG; produced by `spacegraphcats.catlas.catlas`
    * catlas.csv - the catlas for the cDBG; produced by `spacegraphcats.catlas.catlas`
    * commands.log - a partial log of all of the commands
* the `dory_k21_r1_search_oh0/` directory contains the output of a search:
    * results.csv - summary results for the queries (containment, similarity, etc.)
    * dory-head.fa.cdbg_ids.txt.gz - cDBG node IDs (unitig IDs) matching query
    * dory-head.fa.contigs.sig - the sourmash signature of the entire match in the cDBG
    * dory-head.fa.frontier.txt.gz - catlas node IDs matching query
    * dory-head.fa.response.txt - response curve showing how much overhead is gained for each node
    * command.txt - a partial log of the commands run

