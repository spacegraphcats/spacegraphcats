#! /usr/bin/env python
import os
import sys
import argparse
from sourmash_lib import MinHash
from collections import defaultdict
from spacegraphcats.catlas import CAtlas
import time
import sourmash_lib
from sourmash_lib.sbt import SBT, GraphFactory
from sourmash_lib.sbtmh import search_minhashes, SigLeaf
from sourmash_lib import signature
import screed
import pickle


class MinHashFactory(object):
    def __init__(self, **params):
        self.params = params

    def __call__(self):
        return MinHash(**self.params)


def make_contig_minhashes(contigfile, factory):
    "Make the minhashes for each contig in the contigfile."

    d = {}
    total_bp = 0
    watermark = 1e7
    for record in screed.open(contigfile):
        if total_bp >= watermark:
            print('... {:5.0e} bp thru contigs'.format(int(watermark)),
                  file=sys.stderr)
            watermark += 1e7

        mh = factory()
        mh.add_sequence(record.sequence)
        d[int(record.name)] = mh

        total_bp += len(record.sequence)

    return d


def load_layer1_to_cdbg(catlas_file, domfile):
    "Load the mapping between first layer catlas and the original DBG nodes."

    domset = {}
    for line in open(domfile, 'rt'):
        dom_node, *beneath = line.strip().split(' ')

        dom_node = int(dom_node)
        beneath = map(int, beneath)

        domset[dom_node] = set(beneath)

    layer1_to_cdbg = defaultdict(set)

    catlas_to_cdbg = {}
    for line in open(catlas_file, 'rt'):
        catlas_node, cdbg_node, level, beneath = line.strip().split(',')
        if int(level) != 0:
            continue

        catlas_node = int(catlas_node)
        cdbg_node = int(cdbg_node)
        catlas_to_cdbg[catlas_node] = domset[cdbg_node]

    for line in open(catlas_file, 'rt'):
        catlas_node, cdbg_node, level, beneath = line.strip().split(',')
        if int(level) != 1:
            continue

        catlas_node = int(catlas_node)
        beneath = beneath.split(' ')
        beneath = list(map(int, beneath))

        for catnode in beneath:
            dbgnodes = catlas_to_cdbg[catnode]
            layer1_to_cdbg[catlas_node].update(dbgnodes)

    return layer1_to_cdbg


def build_dag(catlas_file, leaf_minhashes, factory):
    "Build MinHashes for all the internal nodes of the catlas DAG."

    # create a list of all the nodes, sorted by level (increasing)
    x = []
    for line in open(catlas_file, 'rt'):
        catlas_node, cdbg_node, level, beneath = line.strip().split(',')
        if int(level) <= 1:
            continue

        beneath = beneath.split(' ')
        beneath = list(map(int, beneath))

        x.append((int(level), int(catlas_node), beneath))

    # walk through, building the merged minhashes (which we can do in a
    # single pass on the sorted list).
    x.sort()
    for (level, catlas_node, beneath) in x:
        merged_mh = factory()

        for subnode in beneath:
            mh = leaf_minhashes[subnode]
            merged_mh.add_many(mh.get_mins())
        leaf_minhashes[catlas_node] = merged_mh


def merge_nodes(child_dict, child_node_list, factory):
    """Merge child nodes into a single minhash."""
    # merge into a single minhash!
    merged_mh = factory()

    for graph_node in child_node_list:
        if graph_node in child_dict:
            mh = child_dict[graph_node]
            merged_mh.merge(mh)

    # add into merged minhashes table.
    return merged_mh


def main():
    p = argparse.ArgumentParser()
    p.add_argument('catlas_prefix', help='catlas prefix')
    p.add_argument('-x', '--bf-size', type=float, default=1e4)
    p.add_argument('--leaves-only', action='store_true')
    p.add_argument('--scaled', default=100.0, type=float)
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('--sbt', action='store_true', help='build SBT for use with sourmash')
    p.add_argument('--sigs', action='store_true', help='save built minhashes for use with sourmash')
    p.add_argument('-o', '--output', default=None)
    p.add_argument('--track-abundance', action='store_true')

    args = p.parse_args()

    ksize = args.ksize
    scaled = args.scaled

    factory = MinHashFactory(n=0, ksize=ksize,
                             max_hash=sourmash_lib.scaled_to_max_hash(scaled),
                             track_abundance=args.track_abundance)
    
    basename = os.path.basename(args.catlas_prefix)
    contigfile = '%s.gxt.contigs' % (basename,)
    contigfile = os.path.join(args.catlas_prefix, contigfile)

    catlas = os.path.join(args.catlas_prefix, basename + '.catlas')
    domfile = os.path.join(args.catlas_prefix, basename + '.domfile')
    
    # make minhashes from node contigs
    print('ksize={} scaled={:.0f}'.format(ksize, scaled))
    print('making contig minhashes...')
    graph_minhashes = make_contig_minhashes(contigfile, factory)
    print('...made {} contig minhashes'.format(len(graph_minhashes)))

    # load mapping between dom nodes and cDBG/graph nodes:
    layer1_to_cdbg = load_layer1_to_cdbg(catlas, domfile)
    print('loaded {} layer 1 catlas nodes'.format(len(layer1_to_cdbg)))
    x = set()
    for v in layer1_to_cdbg.values():
        x.update(v)
    print('...corresponding to {} cDBG nodes.'.format(len(x)))

    # create minhashes for catlas leaf nodes.
    leaf_minhashes = {}
    for n, (catlas_node, cdbg_nodes) in enumerate(layer1_to_cdbg.items()):
        if n and n % 1000 == 0:
            print('... built {} leaf node MinHashes...'.format(n),
                  file=sys.stderr)
        mh = merge_nodes(graph_minhashes, cdbg_nodes, factory)
        leaf_minhashes[catlas_node] = mh
    print('created {} leaf node MinHashes via merging'.format(n + 1))
    print('')

    # build minhashes for entire catlas, or just the leaves (dom nodes)?
    if not args.leaves_only:
        build_dag(catlas, leaf_minhashes, factory)

    path = os.path.basename(args.catlas_prefix) + '.minhashes'
    path = os.path.join(args.catlas_prefix, path)

    try:
        os.mkdir(path)
    except FileExistsError:
        pass

    print('saving individual minhashes in {}/node*.pickle'.format(path))
    for node_id, mh in leaf_minhashes.items():
        name = 'node{}.pickle'.format(node_id)
        name = os.path.join(path, name)
        with open(name, 'wb') as fp:
            pickle.dump(mh, fp)

    if args.sbt or args.sigs:
        print('')
        print('building signatures for use with sourmash')

        sigs = []
        for node_id, mh in leaf_minhashes.items():
            ss = signature.SourmashSignature('', mh,
                                             name='node{}'.format(node_id))
            sigs.append(ss)

        # shall we output an SBT or just a file full of signatures?
        if args.sigs:
            # just sigs - decide on output name
            if args.output:
                signame = args.output
            else:
                signame = os.path.basename(args.catlas_prefix) + '.sig'
                signame = os.path.join(args.catlas_prefix, signame)

            print('saving sigs to "{}"'.format(signame))

            with open(signame, 'wt') as fp:
                signature.save_signatures(sigs, fp)

        if args.sbt:
            # build an SBT!
            factory = GraphFactory(1, args.bf_size, 4)
            tree = SBT(factory)

            print('')
            print('building Sequence Bloom tree for signatures...')
            for ss in sigs:
                leaf = SigLeaf(ss.md5sum(), ss)
                tree.add_node(leaf)

            print('...done with {} minhashes. saving!'.format(len(leaf_minhashes)))

            if args.output:
                sbt_name = args.output
            else:
                sbt_name = os.path.basename(args.catlas_prefix)
                sbt_name = os.path.join(args.catlas_prefix, sbt_name)
            tree.save(sbt_name)
            print('saved sbt "{}.sbt.json"'.format(sbt_name))

    sys.exit(0)

if __name__ == '__main__':
    main()
