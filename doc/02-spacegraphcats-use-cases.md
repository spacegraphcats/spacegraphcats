# spacegraphcats: use cases and working with the output

## Installing the spacegraphcats software and its dependencies

Please see [Installing spacegraphcats](00-installing-spacegraphcats.md).

## Spacegraphcats use cases

### Metagenome bin completion

As discussed in the [primer on metagenomics and assembly graphs](0a-primer.md), assembly and binning are lossy processes that leave many reads uncharacterized in metagenome analysis. 
Spacegraphcats can be used to recover unassembled and unbinned reads that "belong" to a metagenome bin by using the metagenome bin as a query. 
We recommend the [Atlas](https://metagenome-atlas.readthedocs.io/en/latest/) automated metagenome analysis pipeline to produce the initial bins. 
Then, each bin can be used as a query, producing a query neighborhood that reassociates unassembled and unbinned reads that are graph-adjacent to the query.  
In a metagenome from an oil reservoir, bin completion with spacegraphcats identified strain-variable regions in genes that were present in the bin (e.g. *gyrA* in panel A of the figure below) and identified unbinned genes that augmented the functional content of the metagenome by 13% (panel B of the figure below).

![](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-020-02066-4/MediaObjects/13059_2020_2066_Fig4_HTML.png?as=webp) *Source: [Brown et al. 2020](https://doi.org/10.1186/s13059-020-02066-4)*

Install companion software:

```
conda install sourmash fastp khmer
```

Download data:

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/008/SRR1976948/SRR1976948_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/008/SRR1976948/SRR1976948_2.fastq.gz
```

Adapter trim the data:

```
fastp --in1 SRR1976948_1.fastq.gz \
  --in2 SRR1976948_2.fastq.gz \
  --out1 SRR1976948_1.trim.fastq.gz \
  --out2 SRR1976948_2.trim.fastq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 4 \
  --length_required 31 --correction \
  --json SRR1976948.trim.json \
  --html SRR1976948.trim.html
```

k-mer trim the data:
```
interleave-reads.py SRR1976948_1.trim.fastq.gz SRR1976948_2.trim.fastq.gz | \
        trim-low-abund.py --gzip -C 3 -Z 18 -M 20e9 -V - -o SRR1976948.abundtrim.fq.gz
```

Download the query genome (this is a metagenome bin produced from SRR1976948):
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/508/995/GCA_001508995.1_ASM150899v1/GCA_001508995.1_ASM150899v1_genomic.fna.gz
```

Configuration file, saved as `conf1.yml`.
Note that configuration files may contain multiple input sequences from which to build the cDBG and catlas, or multiple search fasta files. 

```
catlas_base: 'SRR1976948'
input_sequences:
- SRR1976948.abundtrim.fq.gz
ksize: 31
radius: 1
search:
- GCA_001508995.1_ASM150899v1_genomic.fna.gz
searchquick: GCA_001508995.1_ASM150899v1_genomic.fna.gz
```

run spacegraphcats:

```
python -m spacegraphcats run conf1.yml extract_contigs extract_reads --nolock 
```

compare the output sequences:

```
sourmash compute -k 21,31,51 --scaled 2000 -o GCA_001508995.1_ASM150899v1_genomic.sig GCA_001508995.1_ASM150899v1_genomic.fna.gz
sourmash compare -k 31 --csv comp.csv *sig SRR1976948_k31_r1_search_oh0/*sig
```

Publications using this approach are listed below:

1. Brown, C.T., Moritz, D., O’Brien, M.P. et al. Exploring neighborhoods in large metagenome assembly graphs using spacegraphcats reveals hidden sequence diversity. Genome Biol 21, 164 (2020). https://doi.org/10.1186/s13059-020-02066-4
2. Lumian, J.E., Jungblut, A.D., Dillion, M.L. et al. Metabolic Capacity of the Antarctic Cyanobacterium Phormidium pseudopriestleyi That Sustains Oxygenic Photosynthesis in the Presence of Hydrogen Sulfide. Genes 12, 426 (2021). https://doi.org/10.3390/genes12030426 

### Identifying the context of e.g. horizontally transferred genes

Spacegraphcats queries don't need to be an entire genome, but can be any size (greater than *k*). 
For example, one can query with a gene of interest. 
In the example below, we identify antibiotic resistance genes in time series human stool metagenomes using [GROOT](https://github.com/will-rowe/groot), and then use the identified genes as spacegraphcats queries.
Using this method, we can see how the context of antibiotic resistance genes changes over time.

The dataset used below if a time series from a single individual with Crohn's disease, sequenced by the [Integrative Human Microbiome Project (iHMP)](https://ibdmdb.org/).
The stool microbiome from this individual was sequenced at 12 different timepoints over 37 weeks, and the individual was on antibiotics at various times during the that timeframe.

![](https://i.imgur.com/H7JiYvr.png) *Times of sequencing and antibiotic exposure for individual H4017 from the Integrative Human Microbiome Project.*

We'll determine the sequence context of one antibiotic resistance gene in one sample using one radius. 


Install companion software needed for this analysis:
```
conda install sourmash=4.0.0 fastp=0.20.1 bbmap=38.70 khmer=3.0.0a3 groot=1.1.2 samtools=1.12 bandage=0.8.1 blast=2.11.0
```

Then, download the the sequencing data and perform quality control.
This sample was taken at week 25 when the individual was on antibiotics.

```
mkdir -p inputs/raw
wget -O inputs/raw/HSM67VFJ.tar https://ibdmdb.org/tunnel/static/HMP2/WGS/1818/HSM67VFJ.tar
tar xf inputs/raw/HSM67VFJ.tar -C inputs/raw
```

Adapter trim and (lightly) quality trim with fastp:
```
mkdir -p outputs/fastp
fastp --in1 inputs/raw/HSM67VFJ_R1.fastq.gz \
      --in2 inputs/raw/HSM67VFJ_R2.fastq.gz \
      --out1 outputs/fastp/HSM67VFJ_R1.trim.fq.gz \
      --out2 outputs/fastp/HSM67VFJ_R2.trim.fq.gz \
      --detect_adapter_for_pe \
      --qualified_quality_phred 4 \
      --length_required 31 --correction \
      --json outputs/fastp/HSM67VFJ.json \
      --html outputs/fastp/HSM67VFJ.html
``` 

Remove "host" (e.g. human) DNA from the sequencing reads
You will need the masked human k-mer data available [here](https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit?usp=sharing).
Download and save this file to `inputs/host/inputs/host/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz`.
This step can take up to 64GB of ram.
To reduce this, change `-Xmx64g` to a lower number, e.g. `-Xmx16g`.
```
mkdir -p outputs/bbduk
bbduk.sh -Xmx64g t=3 k=31 \
     in=outputs/fastp/HSM67VFJ_R1.trim.fq.gz \
     in2=outputs/fastp/HSM67VFJ_R2.trim.fq.gz \
     out=outputs/bbduk/HSM67VFJ_R1.nohost.fq.gz \
     out2=outputs/bbduk/HSM67VFJ_R2.nohost.fq.gz \
     outm=outputs/bbduk/HSM67VFJ_R1.human.fq.gz \
     outm2=outputs/bbduk/HSM67VFJ_R2.human.fq.gz \
     ref=inputs/host/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz
```

K-mer trim the reads.
This step can take up to 60GB of ram as written.
To reduce this, change `-M 60e9` to a lower number, e.g. `-M 16e9`. 
```
mkdir -p outputs/abundtrim
interleave-reads.py outputs/bbduk/HSM67VFJ_R1.nohost.fq.gz outputs/bbduk/HSM67VFJ_R2.nohost.fq.gz | \
      trim-low-abund.py --gzip -C 3 -Z 18 -M 60e9 -V - -o outputs/abundtrim/HSM67VFJ.abundtrim.fq.gz
```

After quality control is complete, we'll use [Groot](https://github.com/will-rowe/groot) to identify reads that encode antibiotic resistance genes. 
Groot uses a database from which it identifies genes of interest. 
Download the ARG90 database:

```
groot get -d arg-annot
```

Then, index the database. `-p` specifies how many threads to use during indexing.
```
groot index -m arg-annot.90 -i groot-index -w 100 -p 8
```

Run groot on the sequencing reads:
```
mkdir outputs/groot
groot align -i groot-index -f outputs/abundtrim/HSM67VFJ.abundtrim.fq.gz \
  -p 2 -g outputs/groot/HSM67VFJ_arg90.graph > outputs/groot/HSM67VFJ_arg90.bam
```

And generate a report that summarizes the antibiotic resistance reads identified.
`-c` indicates the amount of coverage for a gene to be reported.
```
groot report --bamFile outputs/groot/HSM67VFJ_arg90.bam -c .9 > outputs/groot/HSM67VFJ_arg90_report.txt
```
The `groot report` function produces the following report:
```
argannot~~~(Bla)cfxA5~~~AY769934:28-993 1452    966     776M1D186M3D
argannot~~~(Bla)cfxA~~~U38243:150-1115  1551    966     886M80D
argannot~~~(Bla)CEP-A-44~~~U05885:556-1458      624     903     613M2D287M1D
argannot~~~(Bla)cfxA4~~~AY769933:1-966  1453    966     963M3D
argannot~~~(Tet)TetW~~~AJ222769:3687-5606       154     1920    44D1851M25D
argannot~~~(Bla)cfxA2~~~AF504910:1-966  1606    966     501M6D456M3D
```

The columns specify `"arg", "read_count", "gene_length", "coverage"` of each antibiotic resistance gene in the sequencing sample.

Using the bam file which recorded which reads map to which antibiotic resistance genes, we can estimate the proportion of reads that mapped to any ARG.

```
mapped=$(samtools view -c -F 4 outputs/groot/HSM67VFJ_arg90.bam)
proportion=$(echo "scale=10 ; $mapped / 10000000" | bc)
printf "HSM67VFJ,$proportion" > outputs/groot/HSM67VFJ_proportion.txt
```

We see that for this sample, a small proportion of reads mapped against the ARG database

```
HSM67VFJ,.0010801000
```

Now that we've identified which antibiotic resistance genes are in this sample, we can select one of them and use it as a spacegraphcats query to determine its context.
We'll start by downloading the FASTA sequences for antibiotic resistance genes in the ARG database and indexing those sequences with samtools.

```
mkdir -p inputs/arg_db
wget -O inputs/arg_db/argannot-args.fna https://github.com/will-rowe/groot/raw/master/db/full-ARG-databases/arg-annot-db/argannot-args.fna
samtools faidx inputs/arg_db/argannot-args.fna
```

Then, extract a single sequence that was identified in the sample.

```
mkdir -p outputs/arg90_matches
samtools faidx inputs/arg_db/argannot-args.fna argannot~~~(Bla)cfxA4~~~AY769933:1-966 > outputs/arg90_matches/cfxA4_AY769933.fna
```

Make a configuration file for spacegraphcats, and save it as `outputs/sgc_conf/HSM67VFJ_r1_conf.yml`.

```
catlas_base: HSM67VFJ
input_sequences:
- outputs/abundtrim/HSM67VFJ.abundtrim.fq.gz
ksize: 31
radius: 1
search:
- outputs/arg90_matches/cfxA4_AY769933.fna
searchquick: outputs/arg90_matches/cfxA4_AY769933.fna
```

Run spacegraphcats to extract the sequence context around the antibiotic resistance gene.
```
python -m spacegraphcats run \
     outputs/sgc_conf/HSM67VFJ_r1_conf.yml \
     extract_contigs extract_reads \
     --nolock --outdir=outputs/sgc_arg_queries_r1 \
     --rerun-incomplete 
```

Lastly, use bcalm and bandage to plot the context extracted by spacegraphcats.

```
mkdir -p outputs/bcalm/HSM67VFJ_r1
bcalm -in outputs/sgc_arg_queries_r1/HSM67VFJ_k31_r1_search_oh0/cfxA4_AY769933.fna.cdbg_ids.reads.gz \
     -out-dir outputs/bcalm/HSM67VFJ_r1 \
     -kmer-size 31 -abundance-min 1 \
     -out outputs/bcalm/HSM67VFJ_r1/cfxA4_AY769933.fna.cdbg_ids.reads.gz
```

Download the script to conver the bcalm file to a gfa file
```
wget https://raw.githubusercontent.com/spacegraphcats/2018-paper-spacegraphcats/master/pipeline-analyses/variant_snakemake/convertToGFA.py
chmod 777 convertToGFA.py
python ./convertToGFA.py outputs/bcalm/HSM67VFJ_r1/cfxA4_AY769933.fna.cdbg_ids.reads.gz.unitigs.fa outputs/bcalm/HSM67VFJ_r1/cfxA4_AY769933.fna.cdbg_ids.reads.gz.unitigs.gfa 31
```

Use the bandage command line cli to plot the results with the antibiotic resistance gene highlighted.

```
mkdir -p outputs/bandage
Bandage image outputs/bcalm/HSM67VFJ_r1/cfxA4_AY769933.fna.cdbg_ids.reads.gz.unitigs.gfa \
     outputs/bandage/HSM67VFJ_r1/cfxA4_AY769933.fna.cdbg_ids.reads.gz.unitigs.png \
     --query outputs/arg90_matches/cfxA4_AY769933.fna
```

A snakemake pipeline encoding this workflow is available [here](https://github.com/taylorreiter/2021-sgc-arg).
Running this workflow on all samples in the times series, we see that the context of the antibiotic resistance gene *cfxA4* changes over time.
Additionallly, using different radius sizes, we can see that the amount of sequencing context we get back increases with increasing radius size.

![](https://i.imgur.com/jM3Lets.png) *Sequence context of antibiotic resistance gene cfxA4 over time and at different radiuses.*

## Alpha functionality

While we have some tried and true use cases for spacegraphcats, we have even more "alpha functionality" -- functionality for which the code and rationale exist, and may have even been run quite a few times, but that hasn't been extensively explored or applied to many real biological problems yet. 
We document this functionality below.  

### Querying by sourmash minhash hash values

[Sourmash](https://sourmash.readthedocs.io/en/latest/) enables [rapid comparisons across large sequencing datasets](https://f1000research.com/articles/8-1006) using scaled minhashing. 
Essentially, sourmash takes a sequence, decomposes it into k-mers, transforms those k-mers into a number via a hash function, and subsamples the numbers. 
This generates a minhash signatures, or a compressed representation, of the original sequencing data, thereby allowing for rapid comparisons even against millions of genomes.
We enabled querying by hash value to allow integration of sourmash and spacegraphcats workflows. 

```
spacegraphcats <conf> hashval_query
```

```
spacegraphcats <conf> extract_reads_for_hashvals
```

The `hashval_ksize` parameter can be different from the k-mer size used to build the cDBG.
However, for now you can only specify one in a config file; make duplicate config files with different `hashval_ksize` values to do queries on multiple ksizes.

In practice, we have not found querying only by hash value extremely useful; sourmash (typically) subsamples to 1/1000th or 1/2000th of the k-mers in a sequence.
Subsampling facilitates rapid and accurate comparisons between sequences, but this resolution provides an incomplete picture of the sequence landscape, missing a lot of important sequence context (e.g. genes, etc.).
Potentially more importantly, while querying with a hash value may return reads in the neighborhood of the k-mer represented by that hash value, those reads may or may not assemble.
If they do not assemble, it becomes really hard to identify what the functional/taxonomic identity of that hash value may be.
We realized this after implementing `hashval_query` and `extract_reads_for_hashvals`, and then implemented multifast queries (see next section).

### Multifasta queries

To circumvent assembly issues that may occur if reads in a neighborhood of a hash value do not assemble, we developed a "multifasta query" that transfers the annotation of a gene query to hash values in that gene's neighborhood.
The underlying goal of "multifasta queries" is to annotate hash values with genes that are in their graph neighborhood.

To do this, spacegraphcats builds an index with a set of FASTA records, e.g. gene sequences from a reference genome, and then queries with a set of "interesting" hash values.
Spacegraphcats will output a CSV file containing hashval-to-record-name multimappings, where both hash values and record names may appear multiple times.

To run a multifasta query, the configuration file should include multifasta-specific fields:

```
multifasta_reference:
- GCF_000021665.1_ASM2166v1_cds_from_genomic.fna

multifasta_scaled: 1000
multifasta_query_sig: data/63-os223.sig
```

Where the `multifasta_reference` contains the gene in a query of interest, `multifasta_scaled` indicates the scaled value for sourmash to use, and `multifasta_query_sig` is a sourmash signature file that contains the hash values that represent k-mers of interest.

Then, run spacegraphcats:

```
python -m spacegraphcats twofoo multifasta_query
```

### graphgrep

`graphgrep` takes a different approach to querying the CDB than is implemented with the catlas. 
The goal of `graphgrep` is to simplify the retrieval of contigs and/or reads based on k-mer queries into the cDBG.
To do this, it implements a `--radius` option that does neighborhood retrieval in a way that is different from catlas neighborhood expansion.
Briefly, `--radius` does a local expansion around matching cDBG nodes without using dominators, which means 1) the catlas does not need to be built, which helps address our current memory issues where in some cases, catlas building is the blocking step; 2) the reads/contigs output here are different from (a superset of) the reads output by our current neighborhood query.
The k-mer index is still required for `graphgrep`.

A sister command, `graphgrep_iter` does the search without having a k-mer index into the cDBG.

For large and/or many queries, we expect this code to be (much) slower than spacegraphcats neighborhood retrieval. 
`graphgrep_iter` will also be quite slow for large graphs.

The read index is still required to output reads instead of contig.

The k-mer index and the reads index are the slowest of the indexing operations performed in a full catlas build, so `graphgrep*` doesn't solve these performance issues.
However, in the specific case where k-mer indexing or catlas building can't be done due to RAM limitations, this may be a useful solution.

`graphgrep*` is currently implemented in pull request #372 and is not a part of the main code base. 

Basic usage:

Output contigs to stdout:

```
python -m spacegraphcats.search.graphgrep twofoo-short twofoo-short_k31 twofoo-short_k31_r1 tests/test-data/63.short.fa.gz
```

Output reads to stdout:

```
python -m spacegraphcats.search.graphgrep -R twofoo-short twofoo-short_k31 twofoo-short_k31_r1 tests/test-data/63.short.fa.gz
```

### Querying a nucleotide catlas with protein sequences

Some biological questions may be better answered with protein queries instead of with nucleotide queries.
Protein queries work by building secondary indices on top of the cDBG unitigs, allowing spacegraphcats to anchor protein k-mers to their source DNA unitigs and reads.
To acheive this, the catlas is built from the cDBG as normal.
Then, a second index is built from the unitigs where the unitigs are translated into protein space and the protein k-mers are connected to the unitigs.
These two indices support protein queries.

Below we illustrate some results obtained from running a protein query.

Using the podar metagenome, the following code performs a protein query with Proteiniclasticum ruminis proteome for the assembly `GCF_000701905.1`.
Then, DNA contigs are extracted the DNA contigs for the results to see what other nucleotides they match against (genome-wise and taxonomically).

The goal of this small use case was to see how much of what is known to be in the metagenome based on DNA matches (to the `GCF_000701905.1` genome as well as the rest of genbank) is recovered by proteome search, as well as to evaluate the "false positives" caused by the much more sensitive proteome search.

```
# do protein search
python scripts/query-unitigs-prot.py podarV-prot/GCF_000701905.1_ruminis.faa.gz podarV_k31/bcalm.unitigs.db --output-prefix podarV_prot_ruminis

# extract direct-match contigs to FASTA file 
python -m spacegraphcats.search.extract_contigs --contigs-db podarV_k31/bcalm.unitigs.db podarV_prot_ruminis.nodes.gz -o podarV_prot_ruminis/podarV_prot_ruminis.nodes.contigs.fa

# take matching unitig IDs, inflate to neighborhood
python -m spacegraphcats.search.extract_neighborhoods_by_cdbg_ids podarV_k31_r1 podarV_prot_ruminis.nodes.gz -o podarV_prot_ruminis/podarV_prot_ruminis.nbhd.nodes.gz

# extract sequences for nbhd to FASTA file
python -m spacegraphcats.search.extract_contigs --contigs-db podarV_k31/bcalm.unitigs.db -o podarV_prot_ruminis/podarV_prot_ruminis.nbhd.contigs.fa podarV_prot_ruminis/podarV_prot_ruminis.nbhd.nodes.gz

# calculate signature with sourmash
cd podarV_prot_ruminis
sourmash compute -k 31 --scaled 1000 *.fa

# search the sequences
sourmash gather *.nodes.*.sig ~/genome-grist/outputs.paper/genbank/SRR606249.x.genbank.matches.sig
sourmash gather *.nbhd.*.sig ~/genome-grist/outputs.paper/genbank/SRR606249.x.genbank.matches.sig
```

From this code, we learn the following insights about the behavior of protein queries.

First, let's look at the gather results with and without neighborhood expansion. 

*gather results without neighborhood expansion*
```
== This is sourmash version 3.5.0. ==
== Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==

select query k=31 automatically.
loaded query: podarV_prot_ruminis.nodes.cont... (k=31, DNA)
loading from /home/ctbrown/genome-grist/outputs.paper/genbank/SRR606249.x.genbanloaded 73 signatures from /home/ctbrown/genome-grist/outputs.paper/genbank/SRR60loaded 73 signatures.


overlap     p_query p_match
---------   ------- -------
1.9 Mbp       29.0%   77.3%    GCA_003514505.1 Proteiniclasticum sp....
222.0 kbp      3.3%    6.1%    GCA_000015865.1 Hungateiclostridium t...
170.0 kbp      2.5%    6.3%    GCA_000172575.2 Enterococcus faecalis...
158.0 kbp      2.4%    5.8%    GCA_000016545.1 Caldicellulosiruptor ...
146.0 kbp      2.2%    6.6%    GCA_000019085.1 Thermoanaerobacter ps...
118.0 kbp      1.8%    4.1%    GCA_000008185.1 Treponema denticola A...
118.0 kbp      1.7%    4.2%    GCA_000022325.1 Caldicellulosiruptor ...
116.0 kbp      1.7%    1.2%    GCA_000013645.1 Paraburkholderia xeno...
107.0 kbp      1.6%    5.8%    GCA_000021645.1 Dictyoglomus turgidum...
0.6 Mbp        1.6%    3.2%    GCA_900115135.1 Proteiniclasticum rum...
93.0 kbp       1.4%    3.8%    GCA_013137915.1 Fusobacterium nucleat...
87.0 kbp       1.3%    1.7%    GCA_000022185.1 Chloroflexus aurantia...
86.0 kbp       1.3%    1.4%    GCA_002959695.1 Bacteroides thetaiota...
85.0 kbp       1.3%    2.2%    GCA_000007985.2 Geobacter sulfurreduc...
81.0 kbp       1.2%    1.1%    GCA_000009705.1 Nostoc sp. PCC 7120 =...
79.0 kbp       1.2%    3.6%    GCA_000006985.1 Chlorobaculum tepidum...
79.0 kbp       1.2%    2.9%    GCA_000020465.1 Chlorobium limicola D...
75.0 kbp       1.1%    2.6%    GCA_000156375.1 Desulfovibrio piger A...
75.0 kbp       1.1%    3.7%    GCA_900637325.1 Wolinella succinogene...
73.0 kbp       1.1%    3.7%    GCA_000021565.1 Persephonella marina ...
68.0 kbp       1.0%    3.7%    GCA_000018945.1 Thermotoga neapolitan...
66.0 kbp       1.0%    2.2%    GCA_000020645.1 Pelodictyon phaeoclat...
66.0 kbp       1.0%    1.2%    GCA_000195675.1 Bordetella bronchisep...
63.0 kbp       0.9%    2.3%    GCA_000020225.1 Akkermansia muciniphi...
63.0 kbp       0.9%    1.9%    GCA_003513425.1 Desulfovibrio sp., AS...
60.0 kbp       0.9%    2.1%    GCA_000015125.1 Chlorobium phaeobacte...
57.0 kbp       0.8%    1.2%    GCA_000011965.2 Ruegeria pomeroyi DSS...
57.0 kbp       0.8%    0.9%    GCA_000018565.1 Herpetosiphon auranti...
57.0 kbp       0.8%    1.4%    GCA_000022565.1 Acidobacterium capsul...
56.0 kbp       0.8%    1.2%    GCA_000019785.1 Leptothrix cholodnii ...
53.0 kbp       0.8%    1.9%    GCA_000009145.1 Nitrosomonas europaea...
53.0 kbp       0.8%    1.1%    GCA_000012825.1 Bacteroides vulgatus ...
53.0 kbp       0.8%    1.6%    GCA_001638825.1 Deinococcus radiodura...
52.0 kbp       0.8%    1.0%    GCA_000017325.1 Shewanella baltica OS...
51.0 kbp       0.8%    2.4%    GCA_003054575.1 Zymomonas mobilis sub...
51.0 kbp       0.8%    2.4%    GCA_003450455.1 Thermus sp., ASM345045v1
51.0 kbp       0.7%    2.4%    GCA_000016085.1 Chlorobium phaeovibri...
found less than 50.0 kbp in common. => exiting

found 37 matches total;
the recovered matches hit 74.4% of the query
```

gather results against the neighborhood

```
== This is sourmash version 3.5.0. ==
== Please cite Brown and Irber (2016), doi:10.21105/joss.00027. ==

select query k=31 automatically.
loaded query: podarV_prot_ruminis.nbhd.conti... (k=31, DNA)
loading from /home/ctbrown/genome-grist/outputs.paper/genbank/SRR606249.x.genbanloaded 73 signatures from /home/ctbrown/genome-grist/outputs.paper/genbank/SRR60loaded 73 signatures.


overlap     p_query p_match
---------   ------- -------
2.0 Mbp       23.0%   81.2%    GCA_003514505.1 Proteiniclasticum sp....
318.0 kbp      3.6%    8.7%    GCA_000015865.1 Hungateiclostridium t...
227.0 kbp      2.6%    8.4%    GCA_000172575.2 Enterococcus faecalis...
209.0 kbp      2.4%    7.6%    GCA_000016545.1 Caldicellulosiruptor ...
197.0 kbp      2.2%    8.9%    GCA_000019085.1 Thermoanaerobacter ps...
191.0 kbp      2.1%    2.0%    GCA_000013645.1 Paraburkholderia xeno...
161.0 kbp      1.8%    5.8%    GCA_000022325.1 Caldicellulosiruptor ...
155.0 kbp      1.7%    5.4%    GCA_000008185.1 Treponema denticola A...
141.0 kbp      1.6%    5.7%    GCA_013137915.1 Fusobacterium nucleat...
138.0 kbp      1.6%    7.5%    GCA_000021645.1 Dictyoglomus turgidum...
132.0 kbp      1.5%    1.8%    GCA_000009705.1 Nostoc sp. PCC 7120 =...
127.0 kbp      1.4%    3.4%    GCA_000007985.2 Geobacter sulfurreduc...
127.0 kbp      1.4%    2.5%    GCA_000022185.1 Chloroflexus aurantia...
121.0 kbp      1.4%    1.9%    GCA_002959695.1 Bacteroides thetaiota...
0.6 Mbp        1.2%    3.3%    GCA_900115135.1 Proteiniclasticum rum...
107.0 kbp      1.2%    3.7%    GCA_000156375.1 Desulfovibrio piger A...
105.0 kbp      1.2%    4.8%    GCA_000006985.1 Chlorobaculum tepidum...
105.0 kbp      1.2%    5.2%    GCA_900637325.1 Wolinella succinogene...
101.0 kbp      1.1%    3.7%    GCA_000020465.1 Chlorobium limicola D...
101.0 kbp      1.1%    5.1%    GCA_000021565.1 Persephonella marina ...
99.0 kbp       1.1%    3.3%    GCA_000020645.1 Pelodictyon phaeoclat...
96.0 kbp       1.1%    1.8%    GCA_000195675.1 Bordetella bronchisep...
95.0 kbp       1.1%    3.3%    GCA_000015125.1 Chlorobium phaeobacte...
92.0 kbp       1.0%    5.0%    GCA_000018945.1 Thermotoga neapolitan...
92.0 kbp       1.0%    2.8%    GCA_001638825.1 Deinococcus radiodura...
91.0 kbp       1.0%    2.8%    GCA_003513425.1 Desulfovibrio sp., AS...
90.0 kbp       1.0%    1.4%    GCA_000018565.1 Herpetosiphon auranti...
89.0 kbp       1.0%    1.8%    GCA_000012825.1 Bacteroides vulgatus ...
89.0 kbp       1.0%    3.3%    GCA_000020225.1 Akkermansia muciniphi...
87.0 kbp       1.0%    1.9%    GCA_000011965.2 Ruegeria pomeroyi DSS...
86.0 kbp       0.9%    1.7%    GCA_000019785.1 Leptothrix cholodnii ...
83.0 kbp       0.9%    2.0%    GCA_000022565.1 Acidobacterium capsul...
80.0 kbp       0.9%    1.5%    GCA_000017325.1 Shewanella baltica OS...
79.0 kbp       0.9%    3.8%    GCA_000016085.1 Chlorobium phaeovibri...
77.0 kbp       0.9%    2.8%    GCA_000009145.1 Nitrosomonas europaea...
77.0 kbp       0.9%    3.6%    GCA_003054575.1 Zymomonas mobilis sub...
70.0 kbp       0.8%    1.3%    GCA_000007345.1 Methanosarcina acetiv...
68.0 kbp       0.8%    1.5%    GCA_000010305.1 Gemmatimonas aurantia...
68.0 kbp       0.8%    3.2%    GCA_003450455.1 Thermus sp., ASM345045v1
59.0 kbp       0.7%    3.7%    GCA_000020785.1 Hydrogenobaculum sp. ...
56.0 kbp       0.6%    1.0%    GCA_000018265.1 Salinispora arenicola...
54.0 kbp       0.6%    0.8%    GCA_000196115.1 Rhodopirellula baltic...
found less than 50.0 kbp in common. => exiting

found 42 matches total;
the recovered matches hit 75.3% of the query
```

**Insights**

+ Protein queries lead to **good genome recovery**. 
    + From the gather results, we see that 77/81% of a Proteiniclasticum genome (top match) was recovered.
    + It's a different genome than the source of the proteome. This is a little unexpected but could make sense - we assume that the proteome search recovered sequences that overlap slightly more sequences from the best-matched assembly than with the query itself.
+ Subsequent matches have increasingly distant taxonomic relationship to *P. ruminis*. This is likely evidence that the protein query picks up real matches that are more distantly related than would be recovered with nucleotide queries.
    + The second match, *Hungateiclostridium*, is in the same order as *P. ruminis* - Clostridiales - but not the same family or genus.
    + Third best match is to Enterococcus faecalis which is only in the same phylum, Firmicutes, as our query.
    + Fourth best match is to Caldicellulosiruptor, same order.
    + Fifth best match is Thermoanaerobacter, same order.
    + This match structure holds true for both node search and neighborhood search, suggesting that the protein search is pretty powerful at finding good matches and the neighborhood inflation is icing on the cake, so to speak.
+ About 75% of the resulting contigs are matched. When compared against a nucleotide neighborhood search for *P. ruminis* that contains novel *P. ruminis* sequence, we see that the containment for these sequences is pretty good, suggesting that the proteome search picks out the majority of the genome found with the regular nucleotide spacegraphcats query. 

```
% sourmash search --containment ruminis-genome.sgc-paper.fa.sig podarV*.sig
...
similarity   match
----------   -----
 85.0%       podarV_prot_ruminis.nbhd.contigs.fa
 81.1%       podarV_prot_ruminis.nodes.contigs.fa
```

Expanding to the protein, the containment drops, suggesting that extra sequences occur in the protein neigborhoods.
```
% sourmash search ruminis-genome.sgc-paper.fa.sig podarV*.sig
...
similarity   match
----------   -----
 31.8%       podarV_prot_ruminis.nodes.contigs.fa
 26.0%       podarV_prot_ruminis.nbhd.contigs.fa
```

There is a lot to be done with understanding how protein searches or, and what all they bring in. 
Please keep this in mind if you choose to use this functionality!

## Working with spacegraphcats results

spacegraphcats is a powerful framework to organize and access unassembled reads in metagenomes, but the output is admittedly unsatisfying. 
While spacegraphcats will give you the reads or k-mers in the neighborhood of your query, often times the reads you're most interested in are the ones that are the most difficult to work with -- the ones that:

1. do not assemble
2. do not match any sequences in databases

We have some experiences with working with these kinds of reads and outline some approaches we have taken to working with them in the past. 
This is still an active area of research that can benefit from the creativity of the metagenomics community!

### Try an amino acid assembler

Sometimes reads do not assemble due to excessive strain variation.
Amino acid assemblers reduce strain variation by translating reads into amino acid sequences prior to assembly.
In particular, this reduces strain variation due to third base pair wobble, where in the third nucleotide in a codon can vary without changing the amino acid which it encodes.

In the past, we have used the [PLASS](https://github.com/soedinglab/plass) amino acid assembler somewhat successfully on spacegraphcats query neighborhoods.
However, PLASS seems to produce many potential amino acid sequences at least in part influenced by combinatorial overlaps in potential protein sequences.
This output can be difficult to wade through.

Amino acid assemblers only produce amino acid sequences. 
To map sequencing reads back to the assembly to estimate coverage depth or number of reads assembled, we have used [paladin](https://github.com/ToniWestbrook/paladin) to map nucleotide reads to the amino acid assembly.

### Use read-level analysis tools

Tools like [mifaser](https://bromberglab.org/project/mifaser/) and [GROOT](https://github.com/will-rowe/groot) perform annotation on reads instead of on assemblies. 
These tools may provide insight into the content of a spacegraphcats query neighborhood.
