# Subgraph classification

Goal: use a MinHash sketch to retrieve matching regions of a graph.

Inputs:

* a graph
* a MinHash sketch

Output:

* a region of the graph. This could be any of

 - a set of catlas nodes
 - a set of dominating nodes (level 0 catlas nodes)
 * a set of nodes in the contracted De Bruijn graph
 
Metrics:

the usual,

* maximize true positives and true negatives;
* minimize false positives and false negatives;

True positives are returned cDBG nodes that contain part of the graph
from which the MinHash was constructed.  True negatives are cDBG nodes
that were not returned and do not contain part of the graph from which
the MinHash was constructed.  False positives are cDBG nodes that were
returned and do not contain part of the graph from which the MinHash
was constructed. False negatives are cDBG nodes that were not returned
but do contain part of the graph from which the MinHash was
constructed.

Estimating TP/FP/FN/TN from nodes in the dominating set is OK, and the
definitions would change to returning nodes that contain cDBG nodes.

There may be acceptable tradeoffs for different downstream
applications, of course!

Assumptions:

* we can assume that the best match is connected, although if not, that's
  also fine.

Confounding factors:

* it may not be possible to get 100% true negative; there are
  situations (strains and pan-genomes) where genomes may partially
  overlap.

## Benchmarking data sets

### Synthetic / simple data sets

`data/tr-cross.fa` and `data/tr-loop.fa` are both constructed from two
linear paths, `data/tr-1.fa` and `data/tr-2.fa`.

### Real data sets with ground truth

`data/15genome.fa.gz` consists of 15 separate microbial genomes.  The
`data/15genome.fa.*.sigdump.txt` are MinHash signatures for each of
the genomes.

### Real data sets with approximate ground truth

To be described.
