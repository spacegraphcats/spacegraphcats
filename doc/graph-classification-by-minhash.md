# Subgraph classification

Goal: use a MinHash sketch to retrieve matching regions of a dataset (graph).

Inputs:

* D: a dataset of genome(s)
* k: a fixed size for k-mers 
* B: the De Bruijn graph for D
* G: the contracted De Bruijn graph (cDBG) with MinHash sketches of size H 
generated from B.
* T: a target subset of D; we probably want to think of this as a contiguous 
region of a genome, which will in turn be a connected subgraph of B.   
* M: a MinHash sketch generated from T 

Output:

* a set of nodes in G. 
	- This could be represented implicitly by a set of Catlas nodes if desired. 
 
Metrics:

We want to:

* maximize true positives (TP) and true negatives (TN);
* minimize false positives (FP) and false negatives (FN);

True positives are returned cDBG nodes that contain part of T -- the region 
of the genome (graph) from which the MinHash was constructed. True negatives 
are cDBG nodes that were not returned and have no k-mers in common with T. 
False positives and negatives are defined analogously. 

Initially, we encapsulate our performance in two metrics: 
* sensitivity: 100 * TP / (TP + FN)
* specificity: 100 * TN / (TN + FP)

NB: There may be acceptable tradeoffs for different downstream
applications, of course!
     

Baseline: 

If the Catlas were perfectly preserving information from the cDBG, it would have 
sensitivity and specificity equal to those of G -- that is, representing information 
loss in going from sets of k-mers to MinHash sketches. For example, a false negative
would occur if the sketch associated with a node x in G did not contain any k-mer in 
M, but the set of nodes of B it represented in the contraction included some that were
associated with T. 


## Initial Implementations

Each of these is a simple (encapsulated) implementation of search, which finds a set of cDBG nodes, 
then calculates the TP/TN/FP/FN and reports the statistics above using "known-true" labels on B (from T). 

NB: There is currently code duplication in these scripts, and we should 
standardize/modularize before doing any real testing.

Terminology: In the Catlas, we will refer to the set of level 0 nodes as 'domination nodes,' and 
define a node's 'shadow' to be all the domination nodes below it in the hierarchy. The shadow
of a set is the union of the shadows of its members. 

### Titus: four basic strategies

The script [search-for-domgraph-nodes.py](https://github.com/spacegraphcats/spacegraphcats/blob/master/search-for-domgraph-nodes.py)
provides four search strategies: bestnode, searchlevel, gathermins,
and gathermins2. These can be specified with '--strategy', and all
four can be run on a small but real data set ('acido') with 'make
bench' in the main directory.  All these strategies are implemented
in `spacegraphcats/catlas_reader.py`.

All four strategies choose one or more catlas nodes; evaluation is done
by finding the level0/dominating set nodes underneath those nodes,
transferring the cDBG labels onto those nodes, and then computing
TP/FP/FN/TN.

The first strategy, 'bestnode', is probably the dumbest possible
search algorithm on the Catlas: it searches all of the Catlas nodes
with a provided MinHash query, and finds the single best matching
node.

The strategy 'searchlevel' searches a given level of the catlas
(specified by `--searchlevel`) and finds all nodes with a MH match
of above 5%.

The strategy 'gathermins' searches a given level of the catlas
(specified by `--searchlevel`), finds all matches above 5%, orders
them by match score, and then greedily selects the subset that
contains non-overlapping components of the MinHash sketch.

The strategy 'gathermins2' searches a given level of the catlas
(specified by `--searchlevel`), finds all matches above 5%, selects
all their **subnodes**, orders them by match score, and then greedily
selects the subset that contains non-overlapping components of the
MinHash sketch.

You can see a comparison of the sensitivity & specificity of all four
approaches on the 8 subchunks of the acidobacterium genome
[here](https://github.com/spacegraphcats/spacegraphcats/blob/master/doc/plot-benchmark-sens-spec.ipynb).

### Blair & Felix: 
The script
[search-with-catlas.py](https://github.com/spacegraphcats/spacegraphcats/blob/master/search-with-catlas.py)
implements a basic frontier-refinement search on the Catlas. 

High-level description: 

* This iterative method searches the Catlas top-down, refining a frontier 
F_i - a set of Catlas nodes that are the current best-match for the query. 
* At each iteration, a "worst" member of the frontier c_i is identified for 
"refinement" -- replacement by a (possibly empty) subset of its children in 
the next frontier. 
* After a node is refined, the new frontier is re-evaluated using a "quality score" 
-- if this has declined, the process terminates, and the prior iteration's frontier 
is used to calculate the output. 

More details: 

* Initialize with F_0 and c_0 equal to the root node of the Catlas. 
* The worst member of the frontier is the highest-level node in F_i. If there's more than one
at this level, we take the node whose MinHash sketch has minimum overlap with M.
* When refining node c_i, we calculate which children to include in the next frontier F_{i+1} by: 
	- computing the union U_i  (not a join -- this will be larger than any one sketch) 
of the MinHash sketches of all of c_i's children in the Catlas (call them X = {x_1, x_2,...}) 
	- take the intersection Q_i of U_i and M - this is what we want to cover
	- Repeat the following until Q_i is empty: 
		* greedily pick x_j in X so its sketch has maximum intersection with Q_i
		* delete elements of sketch(x_j) from Q_i
		* add x_j to F_{i+1} (and remove it from X)
* A frontier's quality score is the product of (1/(average-height(F_i)+1)) and the Jaccard similarity 
of M with the union (not join) of the MinHash sketches of all the frontier's members. The average-height 
is the average Catlas level of the members. If score(F_{i+1} < F_{i}), terminate and return F_i. Otherwise, continue.

As currently implemented, this has a known bias towards higher query coverage 
(sensitivity) at the expense of including too many cDBG nodes (lower specificity). 
Anything in quotations can be tuned/adjusted to improve biological relevance and/or 
performance. Default definitions should be documented in the script comments.

NB: this is an anytime method, meaning stopping at any iteration gives a valid output (the 
shadow of the Catlas nodes currently in the frontier)

NB: experimentation with varying quality scores (including expensive options that look at shadows explicitly
make it clear we need to implement the baseline ASAP to see how good we can possibly do. There appear to be 
phase transitions in similarity/specificity tradeoffs (which can be seen by trading jaccard similarity for max_height
related measures in score). 

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
