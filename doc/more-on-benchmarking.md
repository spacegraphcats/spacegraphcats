# Moar Benchmarking

What's a decent framework for doing benchmarking? How do we record results?
How do we choose parameters? How do we choose domset radius?

## Some basic questions

* how do we know if a node is a good match, or part of a good match?

* if the local complexity of a graph increases (more interconnections,
  for example) how do we have to adjust the search algorithm?

* how quantitative can we be about presence/absence of specific orgs?
  (How quantitative does the dom set let us be, first of all?)

* is @ctb even getting the catlas node/dom node/original graph ID stuff
  right?

## Parameters we can vary

* ksize
* size of minhash when building graph

## Things we could try

Recording the number of k-mers that go into each catlas node.
(This is actually kind of done, I believe?)

Balancing the catlas construction so that, at each level in the DAG,
nodes have approximately the same number of k-mers underneath them.
(This may be confounded by the contraction.)  Perhaps a better idea
would be to allow the search algorithm(s) to take into account the
number of k-mers in the graph beneath each node.

Calculating MinHashes for multiple ksizes (for differential
species/strain resolution a la MetaPalette); storing multiple minhashes
per catlas node.

## Scripts to write

Comparing all nodes at a given level to each other; e.g. how much overlap is
there in the dom set?

## Observations

Domsets nodes don't accurately represent abundance in original graph.


