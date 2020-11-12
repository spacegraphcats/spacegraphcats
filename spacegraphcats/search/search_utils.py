import collections
import csv
import os
import sqlite3

import numpy
from screed.screedRecord import Record
from screed.utils import to_str
from sourmash import MinHash

from spacegraphcats.utils.bgzf.bgzf import BgzfReader
from . import MPHF_KmerIndex


def sqlite_get_max_offset(cursor):
    cursor.execute("SELECT max(sequences.offset) FROM sequences")
    (last_offset,) = next(cursor)

    return last_offset


def sqlite_get_offsets(cursor, cdbg_ids):
    seen_offsets = set()
    seen_labels = set()

    cursor.execute("DROP TABLE IF EXISTS cdbg_query")
    cursor.execute("CREATE TEMPORARY TABLE cdbg_query (cdbg_id INTEGER PRIMARY KEY);")

    for label in cdbg_ids:
        cursor.execute("INSERT INTO cdbg_query (cdbg_id) VALUES (?)", (label,))

    cursor.execute(
        "SELECT DISTINCT sequences.offset,sequences.cdbg_id FROM sequences WHERE cdbg_id in (SELECT cdbg_id FROM cdbg_query) ORDER BY offset"
    )

    for n, (offset, label) in enumerate(cursor):
        if offset not in seen_offsets:
            yield offset
        seen_offsets.add(offset)
        seen_labels.add(label)

    seen_labels -= cdbg_ids
    assert not seen_labels  # should have gotten ALL the labels


class GrabBGZF_Random(object):
    def __init__(self, filename):
        self.reader = BgzfReader(filename, "rt")

        ch = self.reader.read(1)

        if ch == ">":
            iter_fn = my_fasta_iter
        elif ch == "@":
            iter_fn = my_fastq_iter
        else:
            raise Exception("unknown start chr {}".format(ch))

        self.iter_fn = iter_fn

    def get_sequence_at(self, pos):
        self.reader.seek(pos)
        record = next(self.iter_fn(self.reader))
        return record


def iterate_bgzf(reader):
    ch = reader.read(1)

    if ch == ">":
        iter_fn = my_fasta_iter
    elif ch == "@":
        iter_fn = my_fastq_iter
    else:
        raise Exception("unknown start chr {}".format(ch))

    reader.seek(0)

    for record, pos in iter_fn(reader):
        yield record, pos


def my_fasta_iter(handle, parse_description=False, line=None):
    """
    Iterator over the given FASTA file handle, returning records. handle
    is a handle to a file opened for reading
    """
    last_start = handle.tell()
    if line is None:
        line = handle.readline()

    while line:
        data = {}

        line = to_str(line.strip())
        if not line.startswith(">"):
            raise IOError(
                "Bad FASTA format: no '>' at beginning of line: {}".format(line)
            )

        if parse_description:  # Try to grab the name and optional description
            try:
                data["name"], data["description"] = line[1:].split(" ", 1)
            except ValueError:  # No optional description
                data["name"] = line[1:]
                data["description"] = ""
        else:
            data["name"] = line[1:]
            data["description"] = ""

        data["name"] = data["name"].strip()
        data["description"] = data["description"].strip()

        # Collect sequence lines into a list
        sequenceList = []
        pos = handle.tell()
        line = to_str(handle.readline())
        while line and not line.startswith(">"):
            sequenceList.append(line.strip())
            pos = handle.tell()
            line = to_str(handle.readline())

        data["sequence"] = "".join(sequenceList)
        yield Record(**data), last_start
        last_start = pos


def my_fastq_iter(handle, line=None, parse_description=False):
    """
    Iterator over the given FASTQ file handle returning records. handle
    is a handle to a file opened for reading

    CTB: this relies on each FASTQ record being exactly 4 lines.
    """
    while 1:
        pos = handle.tell()

        line = handle.readline()
        if not line:
            return
        assert line.startswith("@"), line
        name = to_str(line.strip())[1:]

        line = handle.readline()
        sequence = to_str(line.strip())

        line = handle.readline()
        plus = to_str(line.strip())
        assert plus == "+"

        line = handle.readline()
        quality = to_str(line.strip())

        yield Record(name, sequence, quality=quality), pos


def get_reads_by_cdbg(sqlite_filename, reads_filename, cdbg_ids):
    """
    Given a list of cDBG IDs, retrieve the actual sequences corresponding
    to them by using offsets into a BGZF file.
    """
    # connect to sqlite db
    db = sqlite3.connect(sqlite_filename)
    cursor = db.cursor()

    # open readsfile for random access
    reads_grabber = GrabBGZF_Random(reads_filename)

    ## get last offset in file as measure of progress
    last_offset = sqlite_get_max_offset(cursor)

    # pull out the offsets of all sequences with matches in cdbg_ids.
    for offset in sqlite_get_offsets(cursor, cdbg_ids):
        offset_f = offset / last_offset

        record, xx = reads_grabber.get_sequence_at(offset)
        assert xx == offset

        yield record, offset_f


def get_contigs_by_cdbg_sqlite(db, cdbg_ids):
    """
    Given a list of cDBG IDs, retrieve the actual contig sequences
    corresponding to them from a sqlite database created by
    sort_bcalm_unitigs.
    """
    cursor = db.cursor()

    for cdbg_id in cdbg_ids:
        cdbg_id = int(cdbg_id)
        cursor.execute("SELECT sequence FROM sequences WHERE id=?", (cdbg_id,))

        results = cursor.fetchall()
        assert len(results) == 1
        (seq,) = results[0]
        yield Record(str(cdbg_id), seq)


def contigs_iter_sqlite(contigs_db):
    """
    Yield all the sequences in the contigs database.
    """
    cursor = contigs_db.cursor()

    cursor.execute("SELECT id, sequence FROM sequences")
    for ident, sequence in cursor:
        yield Record(str(ident), sequence)


### MPHF stuff


def load_kmer_index(catlas_prefix):
    "Load kmer index created by search.contigs_by_kmer."
    return MPHF_KmerIndex.from_catlas_directory(catlas_prefix)


def load_cdbg_size_info(catlas_prefix, min_abund=0.0):
    filename = os.path.join(catlas_prefix, "contigs.info.csv")
    with open(filename, "rt") as fp:
        cdbg_kmer_sizes = {}
        cdbg_weighted_kmer_sizes = {}
        r = csv.DictReader(fp)
        for row in r:
            contig_id = int(row["contig_id"])
            n_kmers = int(row["n_kmers"])
            mean_abund = float(row["mean_abund"])
            if not min_abund or mean_abund >= min_abund:
                cdbg_kmer_sizes[contig_id] = n_kmers
                cdbg_weighted_kmer_sizes[contig_id] = mean_abund * n_kmers

    return cdbg_kmer_sizes, cdbg_weighted_kmer_sizes


def decorate_catlas_with_kmer_sizes(
    layer1_to_cdbg, dag, dag_levels, cdbg_kmer_sizes, cdbg_weighted_kmer_sizes
):
    x = []
    for (node_id, level) in dag_levels.items():
        x.append((level, node_id))
    x.sort()

    node_kmer_sizes = {}
    node_weighted_kmer_sizes = {}
    for level, node_id in x:
        if level == 1:  # aggregate across cDBG nodes
            total_kmers = 0
            total_weighted_kmers = 0
            for cdbg_node in layer1_to_cdbg.get(node_id):
                total_kmers += cdbg_kmer_sizes.get(cdbg_node, 0)
                total_weighted_kmers += cdbg_weighted_kmer_sizes.get(cdbg_node, 0)

            node_kmer_sizes[node_id] = total_kmers
            node_weighted_kmer_sizes[node_id] = total_weighted_kmers
        else:  # aggregate across children
            sub_size = 0
            sub_weighted_size = 0
            for child_id in dag[node_id]:
                sub_size += node_kmer_sizes[child_id]
                sub_weighted_size += node_weighted_kmer_sizes[child_id]
            node_kmer_sizes[node_id] = sub_size
            node_weighted_kmer_sizes[node_id] = sub_weighted_size

    return node_kmer_sizes, node_weighted_kmer_sizes


def output_response_curve(outname, match_counts, kmer_idx, layer1_to_cdbg):
    curve = []

    # track total containment
    total = 0

    # walk over all layer1 nodes
    for node_id, cdbg_nodes in sorted(layer1_to_cdbg.items()):
        n_matches = 0
        n_kmers = 0

        # aggregate counts across cDBG nodes under this layer1 node.
        for cdbg_node in cdbg_nodes:
            n_matches += match_counts.get(cdbg_node, 0)
            n_kmers += kmer_idx.get_cdbg_size(cdbg_node)

        # do we keep this layer1 node, i.e. does it have positive containment?
        if n_matches:
            n_cont = n_matches
            n_oh = n_kmers - n_matches

            total += n_cont
            curve.append((n_cont, n_oh, node_id))

    # sort by absolute containment
    curve.sort(reverse=True)

    # track total containment etc
    sofar = 0
    total_oh = 0
    total_cont = 0

    # CTB: remove redundant sum_cont fields
    # CTB: switch to CSV output
    # CTB: ask Mike what he wants here :)

    with open(outname, "wt") as fp:
        fp.write(
            "sum_cont relative_cont relative_overhead sum_cont2 sum_oh catlas_id\n"
        )

        # only output ~200 points
        sampling_rate = max(int(len(curve) / 200), 1)

        # do the output thing
        for pos, (n_cont, n_oh, node_id) in enumerate(curve):
            sofar += n_cont
            total_oh += n_oh
            total_cont += n_cont

            # output first and last points, as well as at sampling rate.
            if pos % sampling_rate == 0 or pos == 0 or pos + 1 == len(curve):
                fp.write(
                    "{} {} {} {} {} {}\n".format(
                        sofar,
                        total_cont / total,
                        total_oh / total,
                        n_cont,
                        n_oh,
                        node_id,
                    )
                )
