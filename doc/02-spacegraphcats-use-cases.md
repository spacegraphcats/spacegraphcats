# spacegraphcats: use cases and working with the output

## Installing the spacegraphcats software and its dependencies

Please see [Installing spacegraphcats](00-installing-spacegraphcats.md).

## Spacegraphcats use cases

### Metagenome bin completion

### Query by hashval

### Query by gene

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


