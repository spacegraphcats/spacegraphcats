#! /usr/bin/env python
"""
Given a catlas and a multifasta pickle annotation of cDBG output by 
spacegraphcats.search.index_cdbg_by_multifasta, output annotations
by dominator node.
"""
import sys
import argparse
import time
from collections import defaultdict
import csv
import pickle

import spacegraphcats
from spacegraphcats.utils.logging import notify, error, debug
from spacegraphcats.search import search_utils
from spacegraphcats.search import MPHF_KmerIndex, hash_sequence
from spacegraphcats.search.catlas import CAtlas
import screed


def main():
    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('multifasta_pickle')
    p.add_argument('cdbg_output')
    p.add_argument('dom_output')
    args = p.parse_args()
    
    # load catlas DAG
    catlas = CAtlas(args.catlas_prefix, load_sizefile=True)
    notify('loaded {} nodes from catlas {}', len(catlas), args.catlas_prefix)
    notify('loaded {} layer 1 catlas nodes', len(catlas.layer1_to_cdbg))

    ki_start = time.time()
    kmer_idx = MPHF_KmerIndex.from_catlas_directory(args.catlas_prefix)
    notify('loaded {} k-mers in index ({:.1f}s)',
           len(kmer_idx.mphf_to_kmer), time.time() - ki_start)

    notify(f'loading multifasta pickle from {args.multifasta_pickle}')
    with open(args.multifasta_pickle, 'rb') as fp:
        catlas_prefix, records_to_cdbg, cdbg_to_records = pickle.load(fp)

    assert catlas_prefix == args.catlas_prefix, (catlas_prefix, catlas_prefix)

    dom_annots = defaultdict(set)
    for cdbg_id, annots in cdbg_to_records.items():
        dom_id = catlas.cdbg_to_layer1[cdbg_id]
        dom_annots[dom_id].update(annots)

    notify(f'outputting cDBG annots to {args.cdbg_output}')
    with open(args.cdbg_output, 'wt') as fp:
        w = csv.writer(fp)
        w.writerow(['cdbg_id', 'catlas_base', 'filename', 'record_name'])
        for cdbg_id in sorted(catlas.cdbg_sizes):   # get all cDBG IDs
            for (filename, annot) in cdbg_to_records.get(cdbg_id, ()):
                w.writerow([cdbg_id, catlas_prefix, filename, annot])
            
    notify(f'outputting dom node annots to {args.dom_output}')
    with open(args.dom_output, 'wt') as fp:
        w = csv.writer(fp)
        w.writerow(['dom_id', 'catlas_base', 'filename', 'record_name'])
        for node_id in sorted(catlas):
            for (filename, annot) in dom_annots.get(node_id, ()):
                w.writerow([node_id, catlas_prefix, filename, annot])
            
    return 0


if __name__ == '__main__':
    sys.exit(main())
