#! /usr/bin/env python
"""
Query a catlas with a sequence (read, contig, or genome), and retrieve
cDBG node IDs and MinHash signatures for the matching unitigs in the graph.

XX optionally retrieve/save contig sequences? It wouldn't add much time...
Just disk space.
"""
import argparse
import csv
import gzip
import os
import sys
import time
import sqlite3

import screed
import sourmash
from sourmash import MinHash
import argparse
import os
import sys
import time
import gzip

from .search_utils import get_reads_by_cdbg

from ..utils.logging import notify, error
from . import search_utils
from . import MPHF_KmerIndex, hash_sequence
from .catlas import CAtlas


def main(argv=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("input_node_list_file", help="a cdbg_ids.txt.gz file")
    p.add_argument("-o", "--output-node-list-file", required=True)
    args = p.parse_args(argv)

    # load catlas DAG
    catlas = CAtlas(args.catlas_prefix)
    notify("loaded {} nodes from catlas {}", len(catlas), args.catlas_prefix)
    notify("loaded {} layer 1 catlas nodes", len(catlas.layer1_to_cdbg))
    catlas_name = os.path.basename(args.catlas_prefix.rstrip(("/")))

    # load cDBG nodes
    print(f"reading cDBG nodes from {args.input_node_list_file}")
    with gzip.open(args.input_node_list_file, "rt") as fp:
        cdbg_ids = set([int(x.strip()) for x in fp if x.strip()])
    print(f"...read {len(cdbg_ids)} nodes total.")

    doms = set()
    for cdbg_id in cdbg_ids:
        # dominator for this cDBG node
        doms.add(catlas.cdbg_to_layer1[cdbg_id])

    print(f"promoted {len(cdbg_ids)} cDBG IDs to {len(doms)} dom nodes.")

    dom_shadow = catlas.shadow(doms)
    print(f"...leading to a neighborhood set of {len(dom_shadow)} cDBG nodes.")

    print(f"({len(dom_shadow - cdbg_ids)} of them new.)")

    with gzip.open(args.output_node_list_file, "wt") as fp:
        id_list = "\n".join([str(x) for x in sorted(dom_shadow)])
        print(id_list, file=fp)

    return 0


if __name__ == "__main__":
    sys.exit(main())

