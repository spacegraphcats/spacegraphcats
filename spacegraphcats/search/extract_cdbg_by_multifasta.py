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
    p.add_argument('--multi-idx')
    p.add_argument('--output')
    args = p.parse_args(argv)

    assert args.multi_idx
    assert args.output

    with open(args.multi_idx, 'rb') as fp:
        catlas_base, records_to_cdbg, cdbg_to_records = pickle.load(fp)
        print(f'loaded {len(cdbg_to_records)} cdbg to sequence record mappings for {catlas_base}')

    with open(args.output, 'wt') as fp:
        w = csv.writer(fp)
        w.writerow(["filename", "record_name", "catlas_base", "cdbg_file", "prospective_reads_file"])

        filenum = 0
        for (filename, record_name), cdbg_ids in records_to_cdbg.items():
            if not cdbg_ids:
                continue
            cdbg_file = f'record{filenum}-nbhd.cdbg_ids.txt.gz'
            cdbg_file = os.path.join(os.path.dirname(args.output), cdbg_file)
            
            reads_file = f'record{filenum}-nbhd.reads.fa.gz'
            reads_file = os.path.join(os.path.dirname(args.output), reads_file)
            with gzip.open(cdbg_file, 'wt') as outfp:
                outfp.write("\n".join([ str(i) for i in cdbg_ids ]))

            w.writerow([filename, record_name, catlas_base, cdbg_file, reads_file])

            filenum += 1

    print(f'wrote {filenum+1} files containing cdbg_ids.')

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
