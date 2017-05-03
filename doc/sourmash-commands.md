# Some sourmash commands

## Example: build a MinHash signature using sourmash.

If you want to take a `.fa` file and build a MinHash signature, use
`sourmash compute`:

    sourmash compute ../data/acido-chunk1.fa.gz -k 31 --scaled 1000
    
and sourmash will put the MinHash signature for all the contents of
`acido-chunk1.fa.gz` in `acido-chunk1.fa.gz.sig`.  This `.sig` file
can then be passed in to most sourmash routines as well as
`search.search_catlas_with_minhash`.

## Example: search leaf minhashes only.

The catlas minhashes are built by `search.make_catlas_minhashes`, and by
default are built and output to `name.minhashes/node*.pickle`.  This
is convenient for `search.search_catlas_with_minhash`, but there are some
alternative approaches that let you output sourmash-compatible
files.

If you want to build a fast search tree for all the leaf signatures, do:

    python -m search.make_catlas_minhashes -k 31 --scaled 1000 acido --sbt --leaves-only

You can turn off the default output with `--no-pickles`, which will
speed things up and also not overwrite existing minhashes, and you can
specify a different location/name for the SBT file with `-o newname`.

In order to search these, you can use `sourmash sbt_search`:

    sourmash sbt_search acido/acido.sbt.json acido-chunk1.fa.gz.sig

You should be able to find the guaranteed best match in the tree with `--best-only`:

    sourmash sbt_search acido/acido.sbt.json acido-chunk1.fa.gz.sig --best-only
    
although at the moment this only finds the match with the most shared hashes.

# Collect optimal (?) set of leaf nodes through sbt_gather

You can collect all of the leaf nodes that disjointly match your query like so:

First, build an SBT of only the leaves:

    python -m search.make_catlas_minhashes -k 31 --scaled 1000 acido --sbt --leaves-only

then, run `sbt_gather`:

    sourmash sbt_gather acido/acido.sbt.json acido/acido-chunk1.fa.gz.sig
