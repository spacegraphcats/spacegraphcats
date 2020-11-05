#! /usr/bin/env python
import argparse
import csv
import sys
import pickle

import sourmash


def main(argv):
    """\
    """

    p = argparse.ArgumentParser(description=main.__doc__)
    p.add_argument("--hashvals")
    p.add_argument("--multi-idx")
    p.add_argument("--query-sig")
    p.add_argument("--output")
    p.add_argument("-k", "--ksize", type=int)
    p.add_argument("--scaled", type=int)
    args = p.parse_args(argv)

    assert args.hashvals
    assert args.multi_idx
    assert args.query_sig
    assert args.output

    with open(args.hashvals, "rb") as fp:
        hashval_to_contig_id = pickle.load(fp)
        print(f"loaded {len(hashval_to_contig_id)} hashval to contig mappings")

    with open(args.multi_idx, "rb") as fp:
        catlas_base, records_to_cdbg, cdbg_to_records = pickle.load(fp)
        print(
            f"loaded {len(cdbg_to_records)} cdbg to sequence record mappings for {catlas_base}"
        )

    query_sig = sourmash.load_one_signature(args.query_sig, ksize=args.ksize)
    mh = query_sig.minhash
    mh = mh.downsample(scaled=args.scaled)

    print(
        f"loaded query sig '{str(query_sig)}' with {len(mh)} hashes at scaled={args.scaled}"
    )

    found_hashvals = set()
    found_records = set()
    found_filenames = set()
    n_rows = 0

    with open(args.output, "wt") as fp:
        w = csv.writer(fp)
        w.writerow(["hashval", "catlas_base", "filename", "record_name"])

        for hashval in mh.hashes:
            cdbg_id = hashval_to_contig_id.get(hashval)
            if cdbg_id:
                record_names = cdbg_to_records[cdbg_id]
                for (filename, name) in record_names:
                    w.writerow([hashval, catlas_base, filename, name])
                    found_records.add(name)
                    found_hashvals.add(hashval)
                    found_filenames.add(filename)
                    n_rows += 1

    print(
        f"wrote {n_rows} rows, with {len(found_records)} distinct records and {len(found_hashvals)} distinct hashvals."
    )
    print(f"used records from {len(found_filenames)} files.")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
