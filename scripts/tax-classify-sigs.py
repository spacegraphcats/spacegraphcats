#! /usr/bin/env python
"""
Calculate the "taxonomic purity" of search results from extract_by_query.
"""
import sys
import screed
import argparse
import sourmash
import numpy
from itertools import product
from sourmash._minhash import hash_murmur
from sourmash.lca import lca_utils
from pickle import dump
import pickle
from tax_classify_utils import summarize_taxonomic_purity


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('sigs', nargs='+')
    parser.add_argument('lca_db')
    args = parser.parse_args()

    minhashes = []
    for filename in args.sigs:
        ss = sourmash.load_one_signature(filename)
        minhashes.append(ss.minhash)

    # load the LCA database
    dblist, ksize, scaled = lca_utils.load_databases([args.lca_db], None)
    db = dblist[0]

    # double check scaled requirements
    some_mh = minhashes[0]
    mh_scaled = some_mh.scaled
    if scaled >= mh_scaled:
        print('** warning: many minhashes will go unclassified because LCA database scaled is {}'.format(scaled), file=sys.stderr)
        print('** warning: the minhash scaled is {}'.format(mh_scaled), file=sys.stderr)

    summarize_taxonomic_purity(minhashes, db, verbose=True, filenames=args.sigs)


if __name__ == '__main__':
    main()
