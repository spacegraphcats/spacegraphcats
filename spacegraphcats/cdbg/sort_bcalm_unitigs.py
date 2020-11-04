#! /usr/bin/env python
"""
Sort the bcalm unitigs.fa output (a cDBG) into deterministic order.

Outputs a mapping rather than a new file!

Also saves all of the graph information for later perusal.

Also outputs a sourmash scaled=1000 signature for the input unitigs.
"""
import sys
import argparse
import pickle
import collections

import sourmash
from spacegraphcats.search.search_utils import my_fasta_iter


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("bcalm_unitigs")
    parser.add_argument("mapping_pickle_out")
    parser.add_argument("-k", "--ksize", type=int, default=31)
    args = parser.parse_args(argv)

    unitigs = args.bcalm_unitigs
    ksize = args.ksize

    empty_mh = sourmash.MinHash(n=1, ksize=ksize)

    # read & find minimum hash:
    total_bp = 0
    hashval_to_cdbg = {}
    neighbors = collections.defaultdict(set)
    mean_abunds = {}
    sizes = {}

    # record input k-mers in a minhash
    in_mh = sourmash.MinHash(0, ksize, scaled=1000)

    print(f"loading BCALM unitigs from {unitigs}...")

    unitigs_fp = open(unitigs, "rt")
    for record, offset in my_fasta_iter(unitigs_fp):
        total_bp += len(record.sequence)

        # first, get unitig ID
        name = record.name
        name_split = name.split()
        contig_id = int(name_split[0])

        # second, track the various links
        links = [x for x in name_split[1:] if x.startswith("L:")]
        link_ids = [x.split(":")[2] for x in links]
        link_ids = [int(x) for x in link_ids if int(x) != contig_id]

        neighbors[contig_id].update(link_ids)

        # third, get mean abund
        abund = [x for x in name_split[1:] if x.startswith("km:")]
        assert len(abund) == 1, abund
        abund = abund[0].split(":")
        assert len(abund) == 3
        abund = float(abund[2])

        # fourth, record the abundance sequence and size
        mean_abunds[contig_id] = abund
        sizes[contig_id] = len(record.sequence) - ksize + 1

        # fifth, get the min hash val for this sequence
        mh = empty_mh.copy_and_clear()
        mh.add_sequence(record.sequence)
        assert len(mh) == 1, (len(mh), record.sequence)
        min_hashval = next(iter(mh.hashes))

        assert min_hashval not in hashval_to_cdbg
        hashval_to_cdbg[min_hashval] = (contig_id, offset)

        # sixth, record input k-mers to a bulk signature
        in_mh.add_sequence(record.sequence)

    print(f"...read {len(hashval_to_cdbg)} unitigs, {total_bp:.2e} bp.")

    # check links -- make sure that source is always in its neighbors edges.
    # (this is a check for a recurring bcalm bug that has to do with some
    # kind of threading problem)
    print("validating link structure...")
    fail = False
    for source in neighbors:
        for nbhd in neighbors[source]:
            if source not in neighbors[nbhd]:
                print(f"{source} -> {nbhd}, but not {nbhd} -> {source}")
                fail = True
    print("...done!")

    if fail:
        return -1

    # sort contigs based on min_hashval!
    print("remapping cDBG IDs...")
    hashval_to_cdbg_items = sorted(hashval_to_cdbg.items())

    # remap everything into new coordinate space
    remapping = {}

    # remap sequences
    new_key = 0
    for hashval, (old_key, offset) in hashval_to_cdbg_items:
        remapping[old_key] = new_key
        new_key += 1
    assert len(remapping) == len(hashval_to_cdbg_items)

    # remap other things
    new_neighbors = collections.defaultdict(set)
    for old_key, vv in neighbors.items():
        new_vv = [remapping[v] for v in vv]
        new_neighbors[remapping[old_key]] = set(new_vv)

    new_mean_abunds = {}
    for old_key, value in mean_abunds.items():
        new_mean_abunds[remapping[old_key]] = value

    new_sizes = {}
    for old_key, value in sizes.items():
        new_sizes[remapping[old_key]] = value

    print("...done!")

    mean_abunds = new_mean_abunds
    sizes = new_sizes
    neighbors = new_neighbors

    ## save!
    print(f"saving mappings to '{args.mapping_pickle_out}'")
    offsets = [offset for (_, (_, offset)) in hashval_to_cdbg_items]
    with open(args.mapping_pickle_out, "wb") as fp:
        pickle.dump((ksize, offsets, neighbors, mean_abunds, sizes), fp)

    # output sourmash signature for input contigs
    in_sig = sourmash.SourmashSignature(in_mh, filename=args.bcalm_unitigs)
    sourmash.save_signatures([in_sig], open(args.bcalm_unitigs + ".sig", "wt"))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
