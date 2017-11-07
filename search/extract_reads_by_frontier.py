#! /usr/bin/env python
"""
Do a frontier search and retrieve reads that match to the cDBG nodes
in the frontier.

Uses the output of label_cdbg to do its magic.

See also extract_contigs_by_frontier for a simpler version that just pulls
out linear segments.
"""
import argparse
import os
import sys
import leveldb

from sourmash_lib.sourmash_args import load_query_signature

from .search_utils import get_minhashdb_name, get_reads_by_cdbg
from spacegraphcats.logging import log
from search.frontier_search import (frontier_search, compute_overhead, find_shadow)
from . import search_utils
from .search_utils import (load_dag, load_layer0_to_cdbg, load_minhash)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query_sig', help='query minhash')
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('overhead', help='\% of overhead', type=float)
    p.add_argument('readsfile')
    p.add_argument('labeled_reads_sqlite')
    p.add_argument('output')
    p.add_argument('--no-empty', action='store_true')
    p.add_argument('--purgatory', action='store_true')
    p.add_argument('--fullstats', action='store_true')
    p.add_argument('-k', '--ksize', default='31', type=str,
                   help='k-mer size (default: 31)')
    p.add_argument('--no-remove-empty', action='store_true')
    p.add_argument('--boost-frontier', action='store_true')
    p.add_argument('--scaled', help='downsample query sig, for debug',
                   type=float)

    args = p.parse_args()

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

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
        db_path = get_minhashdb_name(args.catlas_prefix, ksize, args.scaled, 0)
        if not db_path:
            print('** ERROR, minhash DB does not exist for {}'.format(ksize),
                  file=sys.stderr)
            sys.exit(-1)
        print('loading minhashdb:', db_path)
        minhash_db = leveldb.LevelDB(db_path)

        # load query MinHash
        query_sig = load_query_signature(args.query_sig, ksize=ksize,
                                         select_moltype='DNA')
        print('loaded query sig {}'.format(query_sig.name()))

        if args.scaled:
            scaled = int(args.scaled)
            print('downsampling query signature to {}'.format(scaled))
            query_sig.minhash = query_sig.minhash.downsample_scaled(scaled)

        ###

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

            if args.boost_frontier:
                print("boosting frontier shadow")
                boosted_frontier = search_utils.boost_frontier(frontier,
                                                               frontier_mh,
                                                               dag, dag_up,
                                                               minhash_db,
                                                               top_node_id)
                print("...went from {} to {}".format(len(frontier), len(boosted_frontier)))
                frontier = boosted_frontier

        frontier = keep_frontier

        shadow = find_shadow(frontier, dag)
        print("Size of the frontier shadow: {} ({:.1f}%)".format(len(shadow),
                                                             len(shadow) / len(layer0_to_cdbg)* 100))
        if len(shadow) == len(layer0_to_cdbg):
            print('\n*** WARNING: shadow is the entire graph! ***\n')

        query_size = len(query_sig.minhash.get_mins())
        query_bp = query_size * query_sig.minhash.scaled
        print("Size of query minhash: {} (est {:2.1e} bp)".\
                  format(query_size, query_bp))
        minhash_size = len(frontier_mh.get_mins())
        minhash_bp = minhash_size * frontier_mh.scaled
        print("Size of frontier minhash: {} (est {:2.1e} bp); ratio {:.2f}".\
                  format(minhash_size, minhash_bp, minhash_bp / query_bp))

        log(args.catlas_prefix, sys.argv)

        #### extract reads

        for x in shadow:
            cdbg_shadow.update(layer0_to_cdbg.get(x))

    print('loading graph & labels/foo...')
    
    # sql
    dbfilename = args.labeled_reads_sqlite
    assert os.path.exists(dbfilename), 'sqlite file {} does not exist'.format(dbfilename)

    total_bp = 0
    total_seqs = 0
    output_seqs = 0

    # output sequences here:
    outfp = open(args.output, 'wt')

    # track minhash of retrieved reads using original query minhash:
    reads_minhash = query_sig.minhash.copy_and_clear()

    print('running query...')
    reads_iter = get_reads_by_cdbg(dbfilename, args.readsfile, cdbg_shadow)
    for n, (record, offset_f) in enumerate(reads_iter):
        if n % 10000 == 0:
            print('...at n {} ({:.1f}% of file)'.format(n, offset_f * 100),
                  end='\r')

        # track retrieved minhash
        reads_minhash.add_sequence(record.sequence, True)

        # output!
        outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
        total_bp += len(record.sequence)
        total_seqs += 1

    print('')
    print('fetched {} reads, {} bp matching frontier.'.format(total_seqs, total_bp))

    print('query inclusion by retrieved reads: ', query_sig.minhash.contained_by(reads_minhash))
    reads_mh_down = reads_minhash.downsample_scaled(frontier_mh.scaled)
    print('frontier inclusion by retrieved reads: ', frontier_mh.contained_by(reads_mh_down))

    if query_mh.contained_by(frontier_mh) != query_mh.contained_by(reads_mh_down):
        print('*** WARNING: reads containment != frontier containment.',
              file=sys.stderr)

    sys.exit(0)


if __name__ == '__main__':
    main()
