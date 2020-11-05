#! /usr/bin/env python
import sys
import argparse
import time
from collections import defaultdict
import csv

import spacegraphcats
from spacegraphcats.utils.logging import notify, error, debug
from spacegraphcats.search import search_utils
from spacegraphcats.search import MPHF_KmerIndex, hash_sequence
from spacegraphcats.search.catlas import CAtlas
import screed


def main():
    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument(
        "-k", "--ksize", default=31, type=int, help="k-mer size (default: 31)"
    )
    args = p.parse_args()

    # load catlas DAG
    catlas = CAtlas(args.catlas_prefix, load_sizefile=True)
    notify("loaded {} nodes from catlas {}", len(catlas), args.catlas_prefix)
    notify("loaded {} layer 1 catlas nodes", len(catlas.layer1_to_cdbg))

    root = catlas.root
    root_cdbg_nodes = set(catlas.shadow([root]))
    print(f"root cdbg nodes: {len(root_cdbg_nodes)}")

    cdbg_id_set = set()
    for record in screed.open(f"{args.catlas_prefix}/contigs.fa.gz"):
        cdbg_id_set.add(int(record.name))

    print(f"cdbg ID set: {len(cdbg_id_set)}")

    print(f"root - unitigs: {len(root_cdbg_nodes - cdbg_id_set)}")
    print(f"unitigs - root: {len(cdbg_id_set - root_cdbg_nodes)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
