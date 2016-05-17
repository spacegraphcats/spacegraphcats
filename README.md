# spacegraphcats

supercool stuff from DDD Barnraising

## Setup

* Install dependencies with `pip install -r requirements.txt`
* [Set up git lfs](https://git-lfs.github.com/) and initialize it with
  `git lfs install`

## Some big gxt/mxt files

Grab this:

    curl -O http://athyra.idyll.org/~t/transfer/spacegraphcats-extract.tar 

(It's about 220 MB.)

## Running the pipeline

For a pair of .gxt/.mxt files (say `eldritch.gxt` and `eldritch.mxt`) create
a matching project folder `eldritch`. The folder name must match the .gxt/.mxt name.

Run `python3 build-catlas.py /path/to/eldritch r` where *r* determines the coverage
radius. Larger *r* will create smaller atlases at the cost of precision and will
take longer to compute. Intermediate computational steps are cached in the project
director and will speed up subsequence catlas-computations.

### CAtlas file structure

An r-catlas (where *r* is the same parameter as above) consists of the files `eldrictch.catlas.*r*.gxt`
and `eldrictch.catlas.*r*.mxt`. The .gxt file contains the DAG structure with the following fields for nodes:
* `id`: A unique node unrelated to the id in the source graph
* `size`: In a leaf node this size reflects the number of vertices of the original graph assigned to this node.
  For internal nodes this value is the sum of it's childrens sizes.
* `vertex`: If this node corresponds to a vertex from the source graph then this field contains that vertex' id.
			For the root of the DAG this field will say `root`, in all other cases `virtual`.
* `level`: The level of this node in the DAG. Leaves are on level 0 and the level increases towards the root.

The .mxt file contains a minhash set for every node of the DAG, referenced by the `id` field.

## References

[Graph Modelling Language](https://en.wikipedia.org/wiki/Graph_Modelling_Language)

[Graph Text (GXT)](https://github.com/spacegraphcats/spacegraphcats/blob/master/spacegraphcats/parser-examples/README.md)

- [Safe and complete contig assembly via omnitigs](http://arxiv.org/abs/1601.02932) - a good introduction to the theory behind genome assembly.
- [Assembly complexity of prokaryotic genomes using short reads](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-21) - limitations of our source data type.
- [An analysis of the feasibility of short read sequencing](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1278949/) - limitations of our source data type.
