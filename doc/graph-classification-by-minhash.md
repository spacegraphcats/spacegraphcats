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

## Assumptions and confounding factors:

* we can assume that the best match is connected, although if not, that's
  also fine.

* sometimes we may have multiple MinHash signatures to search with,
  but not always. Generally speaking, we will not have a comprehensive
  idea of what (else) might be in the graph.
  
* it may not be possible to get 100% true negative; there are
  situations (strains and pan-genomes) where genomes may partially
  overlap.
  
* because of the previous point, we may know (& use) the size of the
  query sequence/graph from which the query MinHash was constructed,
  but the matching graph may be bigger or smaller.
  
* more generally, the maximum size of the query sequence/graph may be
  between 5,000 and 5 billion k-mers.  (The former is the size of some
  viruses or phage; the latter is the size of some genomes.) For now,
  we can restrict ourselves to the less general case of 500,000 to 5
  million k-mers.  The MinHash sizes of both the query and the catlas
  can be adjusted as needed.

## An initial implementation

The script
[search-for-domgraph-nodes.py](https://github.com/spacegraphcats/spacegraphcats/blob/master/search-for-domgraph-nodes.py)
is a simple implementation of search that calculates the TP/TN/FP/FN.

Briefly, `search-for-domgraph-nodes.py`:

* searches all of the catlas nodes with a provided MinHash query, and
  finds the best matching node;

* for the best matching node, the script then finds the level 0
  nodes (the domination graph nodes) underneath it;

* these level0 nodes are then connected to the original graph, and
  the recovered nodes are compared with "known-true" labels on the
  original graph.

This is probably the dumbest possible search algorithm on the catlas,
but it's a start :).

## Benchmarking data sets

### Synthetic / simple data sets

`data/tr-cross.fa` and `data/tr-loop.fa` are both constructed from two
linear paths, `data/tr-1.fa` and `data/tr-2.fa`.

Benchmarking results are [here](https://github.com/spacegraphcats/spacegraphcats/blob/master/doc/benchmark-tr-cross.ipynb).

### Real data sets with ground truth

`data/acido.fa.gz` is an Acidobacterium genome downloaded from NCBI.
`data/acido-chunk?.fa.gz` is that genome divided into 8.  The associated
`*.sig.dump.txt` files are the signatures for each cunk.

Benchmarking results for retrieving the first chunk of 8 are
[here](https://github.com/spacegraphcats/spacegraphcats/blob/master/doc/benchmark-acido-chunks.ipynb).

`data/15genome.fa.gz` consists of 15 separate microbial genomes.  The
`data/15genome.fa.*.sigdump.txt` are MinHash signatures for each of
the genomes.

### Real data sets with approximate ground truth

To be described.
