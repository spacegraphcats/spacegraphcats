#! /usr/bin/env python
"""
Convert a bcalm unitigs.fa output (a cDBG) into spacegraphcats files.

Outputs a GXT file (containing an undirected graph), a BGZF file
containing the sequences, and a .info.csv file containing
the BGZF offset, mean abundance, and length of each contig.

Also outputs a sourmash scaled=1000 signature for the output contigs.

Takes as input the preprocessed info coming from sort_bcalm_unitigs.
"""
import sys
import argparse
from spacegraphcats.utils.bgzf import bgzf
from spacegraphcats.search.search_utils import my_fasta_iter
import logging
import sourmash
import os.path
import pickle
import sqlite3


def end_match(s, t, k, direction="sp"):
    """
    Compares the first or last k-1 bases of strings s,t for a match.
    The direction argument is a two character string whose characters are 'p'
    or 's', where the first character indicates whether to use the prefix or
    suffix of s for comparison and the second character indicates likewise for
    t.
    Returns a boolean 2-tuple whose first entry indicates whether there is a
    match and whose second entry indicates whether the reverse complement was
    necessary for a match.  If there is no match, the second entry will always
    be False.
    """
    r = reverse_complement(t)
    if direction not in ["pp", "ss", "ps", "sp"]:
        raise ValueError("Valid directions are 'pp', 'ss', 'ps', 'sp'")
    if direction[0] == "p":
        s_end = s[: k - 1]
    else:
        s_end = s[1 - k :]
    if direction[1] == "p":
        t_end = t[: k - 1]
        r_end = r[: k - 1]
    else:
        t_end = t[1 - k :]
        r_end = r[1 - k :]
    if s_end == t_end:
        return (True, False)
    elif s_end == r_end:
        return (True, True)
    else:
        return (False, False)


def is_directed_path(x, sequences, neighbors, k):
    assert len(neighbors[x]) == 2, ValueError(
        "is_directed_path requires a degree 2 vertex"
    )
    u, v = neighbors[x]
    seq_u = sequences[u]
    seq_v = sequences[v]
    if end_match(seq_u, seq_v, k, "pp") or end_match(seq_u, seq_v, k, "ss"):
        return False
    else:
        return True


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[x] for x in reversed(seq))


def contract_neighbor(x, u, neighbors, sequences, mean_abunds, sizes, k):
    seq_x = sequences[x]
    seq_u = sequences[u]
    r = reverse_complement(seq_u)
    logging.debug("neighbor {}".format(u))
    logging.debug("seq_u: {}, {}".format(seq_u[: k - 1], seq_u[1 - k :]))
    logging.debug("   rc: {}, {}".format(r[: k - 1], r[1 - k :]))
    # check which direction the match occurs
    match, rc = end_match(seq_x, seq_u, k, "sp")
    if match:
        if rc:
            sequences[x] += reverse_complement(seq_u)[k - 1 :]
        else:
            sequences[x] += seq_u[k - 1 :]
    else:
        # need to check whether reverse complement was used
        if seq_x[: k - 1] == seq_u[1 - k :]:
            sequences[x] = seq_u[: 1 - k] + seq_x
        else:
            r = reverse_complement(seq_u)
            assert seq_x[: k - 1] == r[1 - k :]
            sequences[x] = r[: 1 - k] + seq_x

    total_abunds = mean_abunds[x] * sizes[x] + mean_abunds[u] * sizes[u]
    sizes[x] += sizes[u]
    mean_abunds[x] = total_abunds / sizes[x]
    # remove v from the graph by making u's other neighbor the
    # neighbor of x
    # there is no neighbor if u has degree 1
    if len(neighbors[u]) > 1:
        for y in neighbors[u]:
            if y != x:
                break
        else:
            msg = "{} doesn't have {} as a neighbor".format(u, x)
            raise ValueError(msg)
        neighbors[y].remove(u)
        neighbors[y].add(x)
        neighbors[x].add(y)
    else:
        y = None  # for debug purposes
    neighbors[x].remove(u)
    neighbors[u] = set()
    logging.debug("removed {}, replacing it with {}, {}".format(u, x, y))

    del sequences[u]


def contract_degree_two(
    non_pendants, neighbors, sequences, mean_abunds, sizes, k, removed_nodes
):
    deg_2 = list()
    for v, N in sorted(neighbors.items()):  # do we need sorted here!?
        if v in non_pendants or len(N) == 0:
            continue
        u = list(N)[0]
        neighbors[u].remove(v)
        N.remove(u)
        if len(neighbors[u]) == 2:
            deg_2.append(u)
    logging.debug("len(deg_2): {}".format(len(deg_2)))
    for x in deg_2:
        if len(neighbors[x]) != 2:
            continue
        u, v = list(neighbors[x])
        seq_x = sequences[x]
        seq_u = sequences[u]
        seq_v = sequences[v]
        # if uxv doesn't form a directed path, we can't do anything
        if end_match(seq_u, seq_v, k, "pp")[0] or end_match(seq_u, seq_v, k, "ss")[0]:
            continue
        logging.debug("analyzing {}".format(x))
        logging.debug("seq_x: {}, {}".format(seq_x[: k - 1], seq_x[1 - k :]))

        did_delete = False
        # can only delete u or v if they have low degree and their
        # neighbors have a directed path
        u_deg = len(neighbors[u])
        if u_deg == 1 or (u_deg == 2 and is_directed_path(u, sequences, neighbors, k)):
            contract_neighbor(x, u, neighbors, sequences, mean_abunds, sizes, k)
            non_pendants.remove(u)
            removed_nodes.add(u)
            did_delete = True
        v_deg = len(neighbors[v])
        if v_deg == 1 or (v_deg == 2 and is_directed_path(v, sequences, neighbors, k)):
            contract_neighbor(x, v, neighbors, sequences, mean_abunds, sizes, k)
            non_pendants.remove(v)
            removed_nodes.add(v)
            did_delete = True

        if did_delete:
            logging.debug("did a removal.")
        else:
            logging.debug("no removal.")


class SqliteAsDict:
    """
    Do direct read access into a FASTA file, mimicking a dictionary.

    Supports in-memory __setitem__, too, that will override future gets.
    """

    def __init__(self, db):
        cursor = db.cursor()
        cursor.execute("PRAGMA synchronous='OFF'")
        cursor.execute("PRAGMA locking_mode=EXCLUSIVE")

        self.db = db
        self.c = cursor

    def __getitem__(self, key):
        self.c.execute("SELECT sequence FROM sequences WHERE id=?", (key,))

        (seq,) = next(iter(self.c))
        return seq

    def __setitem__(self, key, val):
        self.c.execute("UPDATE sequences SET sequence=? WHERE id=?", (val, key))
        assert self.c.rowcount == 1

    def __delitem__(self, key):
        print(f"removing sequence {key}")
        self.c.execute("DELETE FROM sequences WHERE id=?", (key,))

    def close(self):
        self.db.commit()


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("bcalm_unitigs")
    parser.add_argument("sqlite_db")
    parser.add_argument("mapping_pickle")
    parser.add_argument("gxt_out")
    parser.add_argument("contigs_prefix")
    parser.add_argument("-d", "--debug", action="store_true")
    parser.add_argument(
        "-P",
        "--pendants",
        action="store_true",
        help="don't remove low abundance pendants",
    )
    parser.add_argument("-a", "--abundance", nargs="?", type=float, default=1.1)
    parser.add_argument("--randomize", help="randomize cDBG order")  # @CTB
    args = parser.parse_args(argv)

    trim = not args.pendants
    trim_cutoff = args.abundance

    logfile = os.path.join(os.path.dirname(args.gxt_out), "bcalm_to_gxt.log")
    if args.debug:
        logging.basicConfig(filename=logfile, filemode="w", level=logging.DEBUG)
    else:
        logging.basicConfig(filename=logfile, filemode="w", level=logging.WARNING)

    logging.debug("starting bcalm_to_gxt2 run.")

    with open(args.mapping_pickle, "rb") as fp:
        (ksize, neighbors) = pickle.load(fp)

    print(f"Found {len(neighbors)} input unitigs.")

    out_mh = sourmash.MinHash(0, ksize, scaled=1000)

    db = sqlite3.connect(args.sqlite_db)
    sequences = SqliteAsDict(db)

    # build lengths & sizes dictionary

    cursor = db.cursor()
    cursor.execute("SELECT id, abund, LENGTH(sequence) FROM sequences")
    mean_abunds = {}
    sizes = {}
    for idx, abund, length in cursor:
        mean_abunds[idx] = abund
        sizes[idx] = length - ksize + 1

    # if we are removing pendants, we need to relabel the contigs so they are
    # consecutive integers starting from 0.  If not, we create dummy data
    # structures to make the interface the same elsewhere in the data

    all_nodes = set(neighbors.keys())
    if trim:
        print("removing pendants...")
        non_pendants = set(
            v
            for v, N in neighbors.items()
            if len(N) > 1 or mean_abunds[v] > trim_cutoff
        )

        removed_nodes = set()
        contract_degree_two(
            non_pendants, neighbors, sequences, mean_abunds, sizes, ksize, removed_nodes
        )
        all_nodes -= removed_nodes
        print(f"...removed {len(removed_nodes)} nodes.")

    aliases = {x: i for i, x in enumerate(sorted(all_nodes))}
    n = len(aliases)

    # calculate sourmash sig of output & renumber in the sqlite database
    print("calculating sourmash signature of output unitigs")

    c2 = db.cursor()
    cursor.execute("SELECT sequence, offset FROM sequences ORDER BY min_hashval")
    for new_key, (sequence, offset,) in enumerate(cursor):
        c2.execute("UPDATE sequences SET id=? WHERE offset=?", (new_key, offset))
        out_mh.add_sequence(sequence)
        assert c2.rowcount == 1

    db.commit()

    # start the gxt file by writing the number of nodes (unitigs))
    print(f"Outputting graph information to '{args.gxt_out}'...")
    gxtfp = open(args.gxt_out, "wt")
    gxtfp.write("{}\n".format(n))

    # write out all of the links, in 'from to' format.
    n_edges = 0
    for v, N in sorted(neighbors.items()):
        for u in sorted(N):
            gxtfp.write("{} {}\n".format(aliases[v], aliases[u]))
            n_edges += 1

    print("...done! {} vertices, {} edges".format(n, n_edges))

    info_filename = args.contigs_prefix + ".info.csv"
    info_fp = open(info_filename, "wt")

    info_fp.write("contig_id,offset,mean_abund,n_kmers\n")
    for v, i in aliases.items():
        info_fp.write("{},{},{:.3f},{}\n".format(i, 0, mean_abunds[v], sizes[v]))

    # output sourmash signature for output contigs
    out_sig = sourmash.SourmashSignature(out_mh, filename=args.contigs_prefix)
    sourmash.save_signatures([out_sig], open(args.contigs_prefix + ".sig", "wt"))

    sequences.close()

    return 0


if __name__ == "__main__":
    main(sys.argv[1:])
