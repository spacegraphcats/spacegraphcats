#! /usr/bin/env python
"""
Do a frontier search, find the shadow, and extract the contigs in theshadow.

Since every cDBG node in the shadow has a contig associated with it, this
is relatively direct and straightforward and involves no indexing, just a
linear pass through the (usually rather small) contigs file.

See extract_reads_by_frontier.py for extracting the reads in the shadow.
"""
import argparse
import os
import sys
import leveldb

from sourmash_lib.sourmash_args import load_query_signature
import screed

from spacegraphcats.logging import log
from search.frontier_search import (frontier_search, compute_overhead, find_shadow)
from . import search_utils
from .search_utils import (load_dag, load_layer0_to_cdbg, load_minhash)
from .search_utils import get_minhashdb_name
import khmer.utils


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query_sig', help='query minhash')
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('overhead', help='\% of overhead', type=float)
    p.add_argument('output')
    p.add_argument('--no-empty', action='store_true')
    p.add_argument('--purgatory', action='store_true')
    p.add_argument('--fullstats', action='store_true')
    p.add_argument('-k', '--ksize', default='31', type=str,
                   help='list of k-mer sizes (default: 31)')
    p.add_argument('--no-remove-empty', action='store_true')

    args = p.parse_args()

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')
    contigfile = os.path.join(args.catlas_prefix, "contigs.fa.gz")

    # load catlas DAG
    top_node_id, dag, dag_up, dag_levels = load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    # load mapping between dom nodes and cDBG/graph nodes:
    layer0_to_cdbg = load_layer0_to_cdbg(catlas, domfile)
    print('loaded {} layer 0 catlas nodes'.format(len(layer0_to_cdbg)))
    x = set()
    for v in layer0_to_cdbg.values():
        x.update(v)
    print('...corresponding to {} cDBG nodes.'.format(len(x)))

    # minhash search stuff: iterate over ksizes.
    ksizes = [int(k) for k in args.ksize.split(',')]

    cdbg_shadow = set()
    for ksize in ksizes:
        db_path = get_minhashdb_name(args.catlas_prefix, ksize, 0, 0)
        if not db_path:
            print('** ERROR, minhash DB does not exist for {}'.format(ksize),
                  file=sys.stderr)
            sys.exit(-1)
        minhash_db = leveldb.LevelDB(db_path)

        # load query MinHash
        query_sig = load_query_signature(args.query_sig, ksize=ksize,
                                         select_moltype='DNA')
        print('loaded query sig {} for k={}'.format(query_sig.name(), ksize))

        frontier, num_leaves, num_empty, frontier_mh = frontier_search(query_sig, top_node_id, dag, minhash_db, args.overhead, not args.no_empty, args.purgatory)

        top_mh = load_minhash(top_node_id, minhash_db)
        query_mh = query_sig.minhash.downsample_max_hash(top_mh)
        top_mh = top_mh.downsample_max_hash(query_sig.minhash)
        print("Root containment: {}".format(query_mh.contained_by(top_mh)))
        print("Root similarity: {}".format(query_mh.similarity(top_mh)))

        print("Containment of frontier: {}".format(query_mh.contained_by(frontier_mh)))
        print("Similarity of frontier: {}".format(query_mh.similarity(frontier_mh)))
        print("Size of frontier: {} of {} ({:.3}%)".format(len(frontier), len(dag), 100 * len(frontier) / len(dag)))
        print("Overhead of frontier: {}".format(compute_overhead(frontier_mh, query_mh)))
        print("Number of leaves in the frontier: {}".format(num_leaves))
        print("Number of empty catlas nodes in the frontier: {}".format(num_empty))
        print("")

        if not args.no_remove_empty:
            print("removing empty catlas nodes from the frontier...")
            nonempty_frontier = search_utils.remove_empty_catlas_nodes(frontier,
                                                                       minhash_db)
            print("...went from {} to {}".format(len(frontier), len(nonempty_frontier)))
            frontier = nonempty_frontier

        ## now get the catlas level 0 nodes...
        shadow = find_shadow(frontier, dag)

        print("Size of the frontier shadow: {}".format(len(shadow)))
        if len(shadow) == len(layer0_to_cdbg):
            print('\n*** WARNING: shadow is the entire graph! ***\n')

        # diagnostic output from search:

        query_size = len(query_sig.minhash.get_mins())
        query_bp = query_size * query_sig.minhash.scaled
        print("Size of query minhash: {} (est {:2.1e} bp)".\
                  format(query_size, query_bp))
        minhash_size = len(frontier_mh.get_mins())
        minhash_bp = minhash_size * frontier_mh.scaled
        print("Size of frontier minhash: {} (est {:2.1e} bp); ratio {:.2f}".\
                  format(minhash_size, minhash_bp, minhash_bp / query_bp))

        log(args.catlas_prefix, sys.argv)

        # add the cdbg_ids to the shadow.
        for x in shadow:
            cdbg_shadow.update(layer0_to_cdbg.get(x))

    # track various things
    n = 0
    total_bp = 0
    output_bp = 0
    watermark_size = 1e7
    watermark = watermark_size
    total_seqs = 0
    output_seqs = 0

    contigs_minhash = query_mh.copy_and_clear()

    outfp = open(args.output, 'wt')

    for record in screed.open(contigfile):
        n += 1
        if total_bp >= watermark:
            print('... {:5.2e} bp thru contigs; {} read, {} written'.format(int(watermark), total_seqs, output_seqs),
                  file=sys.stderr, end='\r')
            watermark += watermark_size

        total_bp += len(record.sequence)
        total_seqs += 1

        contig_id = int(record.name)

        if contig_id in cdbg_shadow:
            khmer.utils.write_record(record, outfp)
            output_seqs += 1
            output_bp += len(record.sequence)

            contigs_minhash.add_sequence(record.sequence, True)

    print('')
    print('... done! {} read, {} written'.format(total_seqs, output_seqs),
          file=sys.stderr)
    print('{:5.2e} bp written (of {:5.2e} read)'.format(output_bp, total_bp))

    print('(@@ these numbers are only for the last ksize)')
    print('query inclusion by retrieved contigs: ', query_mh.contained_by(contigs_minhash))
    print('frontier inclusion by retrieved reads: ', frontier_mh.contained_by(contigs_minhash))

    sys.exit(0)


if __name__ == '__main__':
    # import cProfile
    # cProfile.run('main()', 'search_stats')

    main()
