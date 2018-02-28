#! /usr/bin/env python
"""
Calculate the "taxonomic purity" of bins built by characterize_catlas_regions.
"""
import sys
import screed
import argparse
import sourmash_lib
import numpy
from itertools import product
from sourmash_lib._minhash import hash_murmur
from sourmash_lib.lca import lca_utils
from pickle import dump
import pickle
import collections
from tax_classify_utils import summarize_taxonomic_purity


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('lca_db')
    args = parser.parse_args()

    # load the minhashes calculated by search.characterize_catlas_regions
    group_ident = pickle.load(open(args.filename + '.node_mh', 'rb'))

    # load the LCA database
    dblist, ksize, scaled = lca_utils.load_databases([args.lca_db], None)
    db = dblist[0]

    # double check scaled requirements
    some_mh = next(iter(group_ident.values()))
    mh_scaled = some_mh.scaled
    if scaled >= mh_scaled:
        print('** warning: many minhashes will go unclassified because LCA database scaled is {}'.format(scaled), file=sys.stderr)
        print('** warning: the minhash scaled is {}'.format(mh_scaled), file=sys.stderr)

    summarize_taxonomic_purity(group_ident.values(), db)


if __name__ == '__main__':
    main()
