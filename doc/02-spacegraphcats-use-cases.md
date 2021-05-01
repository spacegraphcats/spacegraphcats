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

Configuration file:

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

### Query by hashval

### Query by gene

### Multifasta

## Working with spacegraphcats results

spacegraphcats is a powerful framework to organize and access unassembled reads in metagenomes, but the output is admittedly unsatisfying. While spacegraphcats will give you the reads or k-mers associated with your query, often times the reads you're most interested in are the ones that are the most difficult to work with -- the ones that:
1) do not assemble
2) do not match any sequences in databases

We have some experiences with working with these kinds of reads and outline some approaches we have taken to working with them in the past. 
This is still an active area of research that can benefit from the creativity of the spacegraphcats community!

### Try an amino acid assembler

PLASS produces an embarrassment of riches that can be difficult to wade through.

### Use read-level analysis tools

While something did not assemble in your sample, if a similar environment has been sequenced in the past, there's a chance that that thing may have assembled before.


