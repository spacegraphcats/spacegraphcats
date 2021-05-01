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

### Neighborhood of e.g. horizontally transferred genes

Spacegraphcats queries don't need to be an entire genome, but can be any size (greater than *k*). 
For example, one can query with a gene of interest. 
In the example below, we identify antibiotic resistance genes in time series human stool metagenomes using [GROOT](https://github.com/will-rowe/groot), and then use the identified genes as spacegraphcats queries.
Using this method, we can see how the context of antibiotic resistance genes change over time.  
 
 
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
This is still an active area of research that can benefit from the creativity of the spacegraphcats community!

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

