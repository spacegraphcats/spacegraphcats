#! /usr/bin/env python
import argparse
import csv
import gzip
import os
import sys
import time
import pickle
from collections import defaultdict

import screed
import sourmash
from sourmash import MinHash

from ..utils.logging import notify, error, debug
from . import search_utils
from .index import MPHF_KmerIndex
from .catlas import CAtlas


def main(argv):
    """\
    """

    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument('--hashvals')
    p.add_argument('--multi-idx')
    p.add_argument('--query-sig')
    p.add_argument('--output')
    p.add_argument('-k', '--ksize', type=int)
    p.add_argument('--scaled', type=int)
    args = p.parse_args(argv)

    assert args.hashvals
    assert args.multi_idx
    assert args.query_sig
    assert args.output

    with open(args.hashvals, 'rb') as fp:
        hashval_to_contig_id = pickle.load(fp)
        print(f'loaded {len(hashval_to_contig_id)} hashval to contig mappings')

    with open(args.multi_idx, 'rb') as fp:
        records_to_cdbg, cdbg_to_records = pickle.load(fp)
        print(f'loaded {len(cdbg_to_records)} cdbg to sequence record mappings')

    query_sig = sourmash.load_one_signature(args.query_sig,
                                            ksize=args.ksize)
    mh = query_sig.minhash
    mh = mh.downsample_scaled(args.scaled)

    print(f"loaded query sig '{query_sig.name()}' with {len(mh)} hashes at scaled={args.scaled}")

    found_hashvals = set()
    found_records = set()
    n_rows = 0

    with open(args.output, 'wt') as fp:
        w = csv.writer(fp)
        w.writerow(["hashval", "record_name"])

        for hashval in mh.get_mins():
            cdbg_id = hashval_to_contig_id.get(hashval)
            if cdbg_id:
                record_names = cdbg_to_records[cdbg_id]
                for name in record_names:
                    w.writerow([hashval, name])
                    found_records.add(name)
                    found_hashvals.add(hashval)
                    n_rows += 1

    print(f"wrote {n_rows} rows, with {len(found_records)} distinct records and {len(found_hashvals)} distinct hashvals.")

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
