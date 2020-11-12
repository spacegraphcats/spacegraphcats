#! /usr/bin/env python
"""
Retrieve nodes by MinHash hashvals, using the indices created by
spacegraphcats.cdbg.index_cdbg_by_minhash.
Accepts as input catlas prefix, pickled hashval index, file with
list of hashvals, and the output. Also accepts flags -k for
k-mer size (default 31), --scaled for scaled value for sourmash
signature output (default 1000), and -v for verbose.
"""
import argparse
import csv
import gzip
import os
import sys
import time
import pickle
import sqlite3

import sourmash
from sourmash import MinHash

from . import search_utils
from .catlas import CAtlas
from ..utils.logging import notify, error


class QueryOutput:
    def __init__(self, query_hashval, catlas, cdbg_shadow, mh=None):
        self.query_hashval = query_hashval
        self.catlas = catlas
        self.cdbg_shadow = cdbg_shadow
        self.total_bp = 0
        self.total_seq = 0
        self.contigs = []
        self.mh = mh  # if set, track MinHash of contigs

    def __add_sequence(self, sequence):
        self.total_bp += len(sequence)
        self.total_seq += 1

    def retrieve_contigs(self, contigs_db):
        "extract contigs using cDBG shadow."

        # node list.
        # track extracted info
        retrieve_start = time.time()
        # walk through the contigs, retrieving.
        print("extracting contigs...")
        for n, record in enumerate(
            search_utils.get_contigs_by_cdbg_sqlite(contigs_db, self.cdbg_shadow)
        ):
            self.total_seq += 1
            self.total_bp += len(record.sequence)

            if n and n % 10000 == 0:
                offset_f = self.total_seq / len(self.cdbg_shadow)
                notify(
                    "...at n {} ({:.1f}% of shadow)",
                    self.total_seq,
                    offset_f * 100,
                    end="\r",
                )

            self.contigs.append((record.name, record.sequence))
            if self.mh is not None:
                self.mh.add_sequence(record.sequence, force=True)

        # done - got all contigs!
        notify(
            "...fetched {} contigs, {} bp matching hashval domset. " " ({:.1f}s)",
            self.total_seq,
            self.total_bp,
            time.time() - retrieve_start,
        )

    def write(self, csv_writer, csvoutfp, outdir):
        hashval = self.query_hashval
        bp = self.total_bp
        seqs = self.total_seq

        # output to results.csv!
        csv_writer.writerow([hashval, bp, seqs])
        csvoutfp.flush()

        # TR add contigs folder
        # write out cDBG IDs
        q_name = str(hashval)
        cdbg_listname = os.path.basename(q_name) + ".cdbg_ids.txt.gz"
        with gzip.open(os.path.join(outdir, "contigs", cdbg_listname), "wt") as fp:
            fp.write("\n".join([str(x) for x in sorted(self.cdbg_shadow)]))

        # write out contigs
        contigs_outname = os.path.basename(q_name) + ".contigs.fa.gz"
        with gzip.open(os.path.join(outdir, "contigs", contigs_outname), "wt") as fp:
            for name, sequence in self.contigs:
                fp.write(">{}\n{}\n".format(name, sequence))

        # save minhash?
        if self.mh:
            ss = sourmash.SourmashSignature(
                self.mh, name="hashval query:{}".format(q_name)
            )

            sigfile = os.path.join(outdir, "contigs", q_name + ".contigs.sig")
            with open(sigfile, "wt") as fp:
                sourmash.save_signatures([ss], fp)


def execute_query(hashval, catlas, hashval_to_contig_id, mh=None):
    cdbg_node = hashval_to_contig_id.get(hashval, None)

    if not cdbg_node:
        return None

    # find the dominating node in the catlas
    catlas_node_id = catlas.cdbg_to_layer1[cdbg_node]

    # calculate the leaves (this _should_ just be catlas_node_id for now)
    leaves = catlas.leaves([catlas_node_id])
    assert leaves == {catlas_node_id}

    # calculate associated cDBG nodes
    cdbg_shadow = catlas.shadow(leaves)

    return QueryOutput(hashval, catlas, cdbg_shadow, mh=mh)


def main(argv):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("catlas_prefix", help="catlas prefix")
    p.add_argument("mh_index_picklefile", help="pickled hashval index")
    p.add_argument("hashval_list", help="file with list of hashvals")
    p.add_argument("output")
    p.add_argument("--contigs-db", required=True)
    p.add_argument(
        "-k", "--ksize", default=31, type=int, help="k-mer size (default: 31)"
    )
    p.add_argument(
        "--scaled",
        default=1000,
        type=float,
        help="scaled value for contigs minhash output",
    )
    p.add_argument("-v", "--verbose", action="store_true")

    args = p.parse_args(argv)

    # create output directory if it doesn't exist.
    outdir = args.output
    notify("putting output in {}", outdir)
    os.makedirs(os.path.join(outdir, "contigs"), exist_ok=True)

    if not os.path.isdir(outdir):
        error("output '{}' is not a directory and cannot be made", outdir)
        sys.exit(-1)

    # load picklefile
    with open(args.mh_index_picklefile, "rb") as fp:
        hashval_to_contig_id = pickle.load(fp)
    notify(
        "loaded {} hash value -> cdbg_id mappings from {}",
        len(hashval_to_contig_id),
        args.mh_index_picklefile,
    )

    # load list of desired hashvals
    hashvals = [int(x.strip()) for x in open(args.hashval_list, "rt")]
    hashvals = set(hashvals)
    notify("loaded {} search hashvalues from {}", len(hashvals), args.hashval_list)

    if not len(hashvals):
        print("No hash values to search!", file=sys.stderr)
        sys.exit(-1)

    # load catlas DAG
    catlas = CAtlas(args.catlas_prefix)
    notify("loaded {} nodes from catlas {}", len(catlas), args.catlas_prefix)
    notify("loaded {} layer 1 catlas nodes", len(catlas.layer1_to_cdbg))

    # find the contigs filename
    contigs_db = sqlite3.connect(args.contigs_db)

    # get a single ksize & scaled
    ksize = int(args.ksize)
    scaled = int(args.scaled)

    # record command line
    with open(os.path.join(outdir, "command.txt"), "wt") as fp:
        fp.write(str(sys.argv))
        fp.write("\n")

    # output results.csv in the output directory:
    csvoutfp = open(os.path.join(outdir, "hashval_results.csv"), "wt")
    csv_writer = csv.writer(csvoutfp)
    csv_writer.writerow(["hashval", "bp", "contigs"])

    # iterate over each query, do the thing.
    n_found = 0
    for hashval in hashvals:
        notify("----")
        notify("QUERY HASHVAL: {}", hashval)

        mh = MinHash(0, ksize, scaled=scaled)
        result = execute_query(hashval, catlas, hashval_to_contig_id, mh=mh)
        notify("done searching!")
        if not result:
            notify("no result for hashval {}", hashval)
            continue

        result.retrieve_contigs(contigs_db)
        result.write(csv_writer, csvoutfp, outdir)

        assert hashval in mh.hashes

        n_found += 1
    # end main loop!

    notify("----")
    notify(
        "Done! Found {} hashvals of {} in {} with k={}",
        n_found,
        len(hashvals),
        args.catlas_prefix,
        ksize,
    )
    notify("Results are in directory '{}'", outdir)

    return 0


if __name__ == "__main__":
    main(sys.argv[1:])
