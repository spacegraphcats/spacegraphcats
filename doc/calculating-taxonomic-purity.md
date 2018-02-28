# Calculating taxonomic purity of extracted bins

This is documentation on how to use `scripts/tax-classify.py` to estimate the
taxonomic purity of bins of k-mers.

## Getting set up.

First, we'll need an LCA database, which combines signatures with taxonomic information;
the command [sourmash lca index](https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-lca-index) builds such a database.

For this we need two things:

* a bunch of signatures, with names;
* the taxonomic information for those signatures, calculated e.g. by [dib-lab/2018-ncbi-lineages](https://github.com/dib-lab/2018-ncbi-lineages);

For the podar data set, the taxonomic lineages are here: [spacegraphcats/data/podar-lineages.csv](https://github.com/spacegraphcats/spacegraphcats/blob/mphf/data/podar-lineage.csv)

To calculate the signatures for the podar data set, we should use a lower scaled value, which
will give higher resolution to the taxonomic classification:

```
sourmash compute -k 31 --scaled=100 {?,??}.fa -f --name-from-first
```

Then, to calculate an LCA database, do:

```
sourmash lca index podar-lineage.csv podar-ref.scaled100.lca.json \
    {?,??}.fa.sig \
    --scaled 100 -C 3 --split-identifiers
```

This will result in a file `podar-ref.scaled100.lca.json` that you can use with `tax-classify.py`.

## Running tax-classify on catlas region output

First, generate the `.node_mh` file from a catlas with the appropriate scaled value;
it should be no smaller than the one used when computing the signatures & building the LCA database, above.

```
python -m search.characterize_catlas_regions twofoo_k31_r1 twofoo.vec --scaled=100
```

Then:

```
scripts/tax-classify.py twofoo.vec podar-ref.scaled100.lca.json
```