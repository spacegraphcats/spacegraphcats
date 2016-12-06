# CAtlas search - Dec 2016 update

One of the main opportunities in spacegraphcats is to use the catlas
structure to find neighborhoods in the contracted De Bruijn Graph.
One way to do this is to find collections of nodes in the dominating
set that contain good matches to a MinHash query.  The catlas consists
of DAG on top of this dominating set, and we have annotated the catlas
nodes with MinHash sketches composed from the MinHash sketches of the
nodes beneath them.

Our current search strategy revolves around finding
[a frontier](https://github.com/spacegraphcats/spacegraphcats/blob/master/doc/graph-classification-by-minhash.md#blair--felix),
using a top-down search strategy.  This has been slow for large
catlases, however.

## A different way to search the catlas - SBT search

From some other work going on in the lab, we now have a
[fast way to search large collections of MinHash sketches](http://ivory.idyll.org/blog/2016-sourmash-sbt.html)
using
[Sequence Bloom trees (SBT)](http://www.nature.com/nbt/journal/v34/n3/full/nbt.3442.html).
SBTs provide a log2 tree-based search of k-mer collections and they
are well suited to finding MinHash sketches.

SBTs are a complement to the hierarchical composition approach we've
been using for the catlas MinHash sketches - rather than composing
bottom sketches by eliminating hashes, each SBT node consists of the
union of nodes beneath it.

I have used SBTs to index the MinHash sketches from the dominating set
nodes, and then implemented
[a bottom-node search strategy](https://github.com/spacegraphcats/spacegraphcats/blob/sbt_search/search-sbt-of-mxt.py#L185).
This strategy so far is very naive: it finds the collection of
dominating set nodes whose MinHash Sketches contain any overlap at all
with the query signature.  It does not make use of the catlas structure at
all.

The results are impressive, however - with fast search (5-10 seconds
per search of the 15genome data set on my laptop) we achieve
[perfect specificity, with near perfect sensitivity in many cases](https://github.com/TheoryInPractice/cats-in-practice/blob/sbt_search/pipeline/someplots.ipynb).
Sensitivity declines with smaller query sketches (which are less
likely to match dom node MinHash sketches) and smaller dominating set
radii.  The main tradeoff with respect to frontier search is that we
*only* find dominating set nodes, and we find hundreds to thousands of
them.

## Improving SBT search

### Bottom-up catlas exploration

We definitely want to take advantage of the catlas structure so as to
generalize from a sketch search and combine dominating set nodes into
catlas nodes.

A particularly simple approach would be to move up the catlas as long
as the match with the query MinHash did not decrease.  More
complicated would be to tolerate some percentage decrease in the sketch
match as you move up the catlas.

It seems obvious that this will be much faster than a top down search.

## Comments on the approach

* We may not need to pre-compute MinHash Sketches for each node in the
  catlas any more.

* The impact of the choice of domination radius on this approach is
  big: first, `r` impacts the sensitivity of the
  search, because smaller radii encompass fewer relevant hashes; and
  second, smaller `r` results in many more dominating set nodes, which
  then increases the number of SBT nodes.

* The SBT can take up a fair amount of disk space, but it's comparable
  to the catlas mxt file (which may be no longer needed).

* Computing the SBT can be somewhat slow for big collections, but it seems
  to be at least as fast as building the catlas MinHash collection.

* An extension of the SBT approach would be to build an SBT on the
  catlas structure.  I thought a fair bit about starting this way and
  ultimately couldn't decide if the benefits were worth the implementation
  hassle, so went with code I already had working.
