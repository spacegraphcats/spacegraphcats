# spacegraphcats ![Test](https://github.com/spacegraphcats/spacegraphcats/workflows/Test/badge.svg) [![codecov](https://codecov.io/gh/spacegraphcats/spacegraphcats/branch/master/graph/badge.svg)](https://codecov.io/gh/spacegraphcats/spacegraphcats) [![DOI](https://zenodo.org/badge/58208221.svg)](https://zenodo.org/badge/latestdoi/58208221)


Explore large, annoying graphs using hierarchies of dominating sets - because
in space, no one can hear you miao!

This is a collaboration between the
[Theory In Practice](https://github.com/TheoryInPractice/) lab at NC
State and the
[Lab for Data Intensive Biology](https://github.com/dib-lab/) at UC
Davis, generously supported by the Moore Foundation's
[Data Driven Discovery Initiative](https://www.moore.org/initiative-strategy-detail?initiativeId=data-driven-discovery).

![spacegraphcats graph](https://github.com/spacegraphcats/spacegraphcats/raw/master/pics/logo.png)

## Installation and execution quickstart

See [installation instructions](https://github.com/spacegraphcats/spacegraphcats/blob/master/doc/installing-spacegraphcats.md) and [the run guide](https://github.com/spacegraphcats/spacegraphcats/blob/master/doc/running-spacegraphcats.md).

For help or support with this software, please
[file an issue on GitHub](https://github.com/spacegraphcats/spacegraphcats/issues). Thank
you!

### Quickstart

There are two quickstart examples available! Please see
[dory-example](https://github.com/spacegraphcats/spacegraphcats-dory-example)
and
[twofoo-example](https://github.com/spacegraphcats/spacegraphcats-twofoo-example). The
latter example includes
[a snakemake Snakefile](https://snakemake.readthedocs.io/en/stable/).

### Notable dependencies

spacegraphcats uses code from
[BBHash](https://github.com/rizkg/BBHash), a C++ library for building
minimal perfect hash functions (Guillaume Rizk, Antoine Limasset,
Rayan Chikhi; see
[Limasset et al., 2017, arXiv](https://arxiv.org/abs/1702.03154), as
wrapped by [pybbhash](https://github.com/dib-lab/pybbhash).

spacegraphcats also uses functionality from
[khmer](https://github.com/dib-lab/khmer/) and
[sourmash](https://github.com/dib-lab/sourmash).

## Citation information

See the biorxiv preprint [Exploring neighborhoods in large metagenome assembly graphs reveals hidden sequence diversity](https://biorxiv.org/content/early/2018/11/05/462788), Brown et al., 2018, doi: https://doi.org/10.1101/462788.

## Pointers to interesting code

### Interesting algorithms

The `rdomset` code for efficently calculating a dominating set of a graph
at a given radius R is in [spacegraphcats/catlas/rdomset.py](https://github.com/spacegraphcats/spacegraphcats/blob/master/spacegraphcats/catlas/rdomset.py).

The graph denoising code for removing low-abundance pendants from
BCALM cDBGs is in function `contract_degree_two` in
[cdbg/bcalm_to_gxt.py](https://github.com/spacegraphcats/spacegraphcats/blob/master/spacegraphcats/cdbg/bcalm_to_gxt.py).

Part of the `indexPieces` code for indexing cDBG nodes by dominating
nodes is
[cdbg/index_contigs_by_kmer.py](https://github.com/spacegraphcats/spacegraphcats/blob/master/spacegraphcats/cdbg/index_contigs_by_kmer.py). The
remainder is implemented in `search`, below.

The `search` code for extracting query neighborhoods is in
[search/extract_nodes_by_query.py](https://github.com/spacegraphcats/spacegraphcats/blob/master/spacegraphcats/search/extract_nodes_by_query.py);
see especially the call to `kmer_idx.count_cdbg_matches(...)`.

### Interesting library functionality

Code for indexing large FASTQ/FASTA read files by cDBG unitig, and
extracting the reads corresponding to individual unitigs from BGZF
files, is available in
[cdbg/label_cdbg.py](https://github.com/spacegraphcats/spacegraphcats/blob/master/spacegraphcats/cdbg/label_cdbg.py)
and
[search/search_utils.py](https://github.com/spacegraphcats/spacegraphcats/blob/master/spacegraphcats/search/search_utils.py),
`get_reads_by_cdbg`, respectively.
