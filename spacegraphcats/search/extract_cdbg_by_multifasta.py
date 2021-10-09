#! /usr/bin/env python
"""
Extract useful information from a multifasta annotation of the cDBG nodes.

This script takes in a multifasta_index produced by index_cdbg_by_multifasta*
and outputs:
* a CSV file with one multifasta query per line, containing summary information
  about the cDBG nodes to which that query matches.
* individual files for each annotation, listing the cDBG IDs and the name
  of a (potential) output file for where reads will be extracted.
"""
import argparse
import csv
import gzip
import os
import sys
import pickle


def main(argv):
    """\
    """

    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument("--multi-idx", required=True)
    p.add_argument("--info-csv", required=True)
    p.add_argument("--output-cdbg-record", required=True)
    p.add_argument("--output-cdbg-annot", required=True)
    args = p.parse_args(argv)

    # load records <-> cDBG dictionaries
    with open(args.multi_idx, "rb") as fp:
        catlas_base, records_to_cdbg, cdbg_to_records = pickle.load(fp)
        print(
            f"loaded {len(cdbg_to_records)} cdbg to sequence record mappings for {catlas_base}"
        )

    # read in the information about the underlying cDBG nodes
    cdbg_info = {}
    with open(args.info_csv, "rt") as fp:
        r = csv.DictReader(fp)
        for row in r:
            cdbg_id = int(row["contig_id"])
            cdbg_info[cdbg_id] = row
        print(f"loaded {len(cdbg_info)} info records for {catlas_base}")

    # create the per-annotation information file.
    with open(args.output_cdbg_record, "wt") as fp:
        w = csv.writer(fp)
        w.writerow(
            [
                "filename",
                "record_name",
                "catlas_base",
                "cdbg_file",
                "prospective_reads_file",
                "total_kmers",
                "mean_abund",
            ]
        )

        # iterate over all annotations {queryfile, query_record} -> cDBG IDs
        filenum = 0
        for (filename, record_name), cdbg_ids in records_to_cdbg.items():
            if not cdbg_ids:
                continue

            # write out cdbg id list
            cdbg_file = f"record{filenum}-nbhd.cdbg_ids.txt.gz"
            cdbg_file = os.path.join(os.path.dirname(args.output_cdbg_record), cdbg_file)
            with gzip.open(cdbg_file, "wt") as outfp:
                outfp.write("\n".join([str(i) for i in cdbg_ids]))

            # generate *name* of reads.fa.gz file to make
            reads_file = f"record{filenum}-nbhd.reads.fa.gz"
            reads_file = os.path.join(os.path.dirname(args.output_cdbg_record), reads_file)

            # generate info: total k-mers, mean_abund / weighted
            total_kmers = 0
            summed_abund = 0
            for cdbg_id in cdbg_ids:
                info = cdbg_info[cdbg_id]

                n_kmers = int(info["n_kmers"])
                mean_abund = float(info["mean_abund"])
                total_kmers += n_kmers
                summed_abund += n_kmers * mean_abund

            average_abund = summed_abund / total_kmers

            w.writerow(
                [
                    filename,
                    record_name,
                    catlas_base,
                    cdbg_file,
                    reads_file,
                    total_kmers,
                    average_abund,
                ]
            )

            filenum += 1

    print(f"wrote {filenum+1} files containing cdbg_ids.")

    # write out cdbg id to record mapping, one per row
    with open(args.output_cdbg_annot, "wt") as fp:
        w = csv.writer(fp)
        w.writerow(
            [
                "filename",
                "record_name",
                "record_number",
                "catlas_base",
                "cdbg_id",
            ]
        )
    
        filenum = 0
        for (filename, record_name), cdbg_ids in records_to_cdbg.items():
            if not cdbg_ids:
                continue

            record_number = f"record{filenum}"
            for cdbg_id in cdbg_ids:
                cdbg_id = str(cdbg_id)
                w.writerow(
                     [
                        filename,
                        record_name,
                        record_number,
                        catlas_base,
                        cdbg_id,
                    ]
                )

            filenum += 1

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
