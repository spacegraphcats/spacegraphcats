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
conda install sourmash=4.0.0 fastp=0.20.1 bbmap=38.70 khmer=3.0.0a3 groot=1.1.2 samtools=1.12 seqkit=0.16.0 bandage=0.8.1 blast=2.11.0
```

Then, download the the sequencing data and perform quality control.
This sample was taken at week 25 when the individual was on antibiotics.

```
mkdir -p inputs/raw
wget -O inputs/raw/HSM67VFJ.tar https://ibdmdb.org/tunnel/static/HMP2/WGS/1818/HSM67VFJ.tar
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
A snakemake pipeline encoding this workflow is available [here](https://github.com/taylorreiter/2021-sgc-arg).
Running this workflow on all samples in the times series, we see that the context of the antibiotic resistance gene *cfxA4* changes over time.

![](https://i.imgur.com/jM3Lets.png) *Sequence context of antibiotic resistance gene cfxA4 over time and at different radiuses.*

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

### Multifasta

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
