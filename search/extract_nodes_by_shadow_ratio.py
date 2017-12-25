#! /usr/bin/env python
"""
Look for catlas nodes that have few k-mers but many cDBG nodes underneath
them, as a likely sign of strain variation.
"""
import argparse
import os
import sys
import collections
import math
import sourmash_lib

import screed
from . import search_utils
from .search_utils import get_minhashdb_name, get_reads_by_cdbg, load_minhash
from .frontier_search import find_shadow


def get_minhash_size(node_id, minhash_db):
    mh = load_minhash(node_id, minhash_db)
    if mh:
        return len(mh.get_mins())
    return 0


def main(args=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('output')
    p.add_argument('--maxsize', type=float, default=1000)
    p.add_argument('-k', '--ksize', default=31, type=int,
                        help='k-mer size (default: 31)')
    args = p.parse_args(args)

    print('maxsize: {:g}'.format(args.maxsize))

    basename = os.path.basename(args.catlas_prefix)
    catlas = os.path.join(args.catlas_prefix, 'catlas.csv')
    domfile = os.path.join(args.catlas_prefix, 'first_doms.txt')

    # load catlas DAG
    top_node_id, dag, dag_up, dag_levels, cdbg_to_catlas = search_utils.load_dag(catlas)
    print('loaded {} nodes from catlas {}'.format(len(dag), catlas))

    # load mapping between dom nodes and cDBG/graph nodes:
    layer1_to_cdbg = search_utils.load_layer1_to_cdbg(cdbg_to_catlas, domfile)
    print('loaded {} layer 1 catlas nodes'.format(len(layer1_to_cdbg)))

    # load minhash DB
    db_path = get_minhashdb_name(args.catlas_prefix, args.ksize, 0, 0, 43)
    if not db_path:
        print('** ERROR, minhash DB does not exist for {}'.format(args.ksize),
              file=sys.stderr)
        sys.exit(-1)
    minhash_db = search_utils.MinhashSqlDB(db_path)

    # calculate the cDBG shadow sizes for each catlas node.
    x = []
    for (node_id, level) in dag_levels.items():
        x.append((level, node_id))
    x.sort()

    node_shadow_sizes = {}
    for level, node_id in x:
        if level == 1:
            node_shadow_sizes[node_id] = len(layer1_to_cdbg[node_id])
        else:
            sub_size = 0
            for child_id in dag[node_id]:
                sub_size += node_shadow_sizes[child_id]
            node_shadow_sizes[node_id] = sub_size

    # find highest nodes with shadow size less than given max_size
    def find_terminal_nodes(node_id, max_size):
        node_list = set()
        for sub_id in dag[node_id]:
            # shadow size
            size = node_shadow_sizes[sub_id]

            if size < max_size:
                node_list.add(sub_id)
            else:
                children = find_terminal_nodes(sub_id, max_size)
                node_list.update(children)
            
        return node_list

    print('finding terminal nodes for {}.'.format(args.maxsize))

    if 1:
        terminal = set()
        for subnode in dag[top_node_id]:
            mh = load_minhash(subnode, minhash_db)
            if mh:
                terminal.update(find_terminal_nodes(subnode, args.maxsize))

    else:
        terminal = find_terminal_nodes(top_node_id, args.maxsize)
        print('...got {}'.format(len(terminal)))

    # now, go through and collect minhashes for each of them, and
    # estimate total k-mers underneath; filter terminal nodes.
    x = []
    merge_mh = load_minhash(top_node_id, minhash_db).copy_and_clear()
    n_merged = 0
    new_node_set = set()
    for node_id in terminal:
        # get minhash size under this node
        mh = load_minhash(node_id, minhash_db)
        if mh:
            mh_size = len(mh.get_mins()) * mh.scaled

            # retrieve shadow size, calculate mh_size / shadow_size.
            shadow_size = node_shadow_sizes[node_id]
            ratio = math.log(mh_size, 2) - math.log(shadow_size, 2)

            level = dag_levels[node_id]

            # track basic info
            x.append((ratio, node_id, shadow_size, mh_size, level))

            # keep those with larger shadow than k-mers (somewhat arbitrary :)
            # specifically: 10 k-mers per node, or fewer.

            # THIS NEVER HAPPENS for small maxsize?!
            if 2**ratio < 10:
                new_node_set.add(node_id)

                n_merged += 1
                merge_mh.merge(mh)
        else:
            shadow_size = node_shadow_sizes[node_id]
            x.append((0, node_id, shadow_size, 0, level))
            if shadow_size >= 3:
                new_node_set.add(node_id)

    terminal = new_node_set

    print('terminal node stats for maxsize: {:g}'.format(args.maxsize))
    print('n merged:', n_merged)
    print('n tnodes:', len(terminal), '; k-mer cover:', len(merge_mh.get_mins()) * merge_mh.scaled)
    top_mh = load_minhash(top_node_id, minhash_db)
    print('total k-mers:', len(top_mh.get_mins()) * top_mh.scaled)

    x.sort(reverse=True)
    if 1:
        for (k, v, a, b, c) in x:
            print('ratio: {:.3f}'.format(2**k), '/ shadow size:', a, '/ kmers:', b, '/ level:', c)

    if 0:
        print("removing empty catlas nodes from the terminal nodes...")
        nonempty_terminal = search_utils.remove_empty_catlas_nodes(terminal,
                                                                   minhash_db)
        print("...went from {} to {}".format(len(terminal), len(nonempty_terminal)))
        terminal = nonempty_terminal

    # build cDBG shadow ID list.
    cdbg_shadow = set()
    terminal_shadow = find_shadow(terminal, dag)
    for x in terminal_shadow:
        cdbg_shadow.update(layer1_to_cdbg.get(x))

    #### extract reads
    print('extracting contigs')
    contigs = os.path.join(args.catlas_prefix, 'contigs.fa.gz')

    # track results as signature
    contigs_mh = sourmash_lib.MinHash(n=0, ksize=args.ksize, seed=42,
                                      scaled=top_mh.scaled)

    total_bp = 0
    total_seqs = 0

    outfp = open(args.output, 'wt')
    for n, record in enumerate(screed.open(contigs)):
        if n and n % 10000 == 0:
            offset_f = total_seqs / len(cdbg_shadow)
            print('...at n {} ({:.1f}% of shadow)'.format(total_seqs,
                  offset_f * 100),
                  end='\r')

        # contig names == cDBG IDs
        contig_id = int(record.name)
        if contig_id not in cdbg_shadow:
            continue

        outfp.write('>{}\n{}\n'.format(record.name, record.sequence))
        contigs_mh.add_sequence(record.sequence)

        # track retrieved sequences in a minhash
        total_bp += len(record.sequence)
        total_seqs += 1

    # done - got all contigs!
    print('')
    print('fetched {} contigs, {} bp.'.format(total_seqs, total_bp))

    print('wrote contigs to {}'.format(args.output))
    with open(args.output + '.sig', 'wt') as fp:
        ss = sourmash_lib.SourmashSignature(contigs_mh)
        sourmash_lib.save_signatures([ss], fp)


if __name__ == '__main__':
    main()
