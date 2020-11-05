#! /usr/bin/env python
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
    p.add_argument("--multi-idx")
    p.add_argument("--info-csv")
    p.add_argument("--output")
    args = p.parse_args(argv)

    assert args.multi_idx
    assert args.output

    with open(args.multi_idx, "rb") as fp:
        catlas_base, records_to_cdbg, cdbg_to_records = pickle.load(fp)
        print(
            f"loaded {len(cdbg_to_records)} cdbg to sequence record mappings for {catlas_base}"
        )

    cdbg_info = {}
    with open(args.info_csv, "rt") as fp:
        r = csv.DictReader(fp)
        for row in r:
            cdbg_id = int(row["contig_id"])
            cdbg_info[cdbg_id] = row
        print(f"loaded {len(cdbg_info)} info records for {catlas_base}")

    with open(args.output, "wt") as fp:
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

        filenum = 0
        for (filename, record_name), cdbg_ids in records_to_cdbg.items():
            if not cdbg_ids:
                continue

            # write out cdbg id list
            cdbg_file = f"record{filenum}-nbhd.cdbg_ids.txt.gz"
            cdbg_file = os.path.join(os.path.dirname(args.output), cdbg_file)
            with gzip.open(cdbg_file, "wt") as outfp:
                outfp.write("\n".join([str(i) for i in cdbg_ids]))

            # generate *name* of reads.fa.gz file to make
            reads_file = f"record{filenum}-nbhd.reads.fa.gz"
            reads_file = os.path.join(os.path.dirname(args.output), reads_file)

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

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
