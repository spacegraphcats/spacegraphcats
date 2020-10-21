#! /usr/bin/env python
"""
Check recorded md5 hashes against new files; see 'make twofoo-test'.
"""
import hashlib
import sys
import argparse
import os

saved_md5 = {}
saved_md5["twofoo"] = {
    "catlas.csv": "00480517b18203ba261b1a76c7332d8e",
    "cdbg.gxt": "3cf4b4182343e37a9025d9783014fe61",
    "contigs.fa.gz.sig": "517134d21bce89d7879cfdfcf46e5b8f",
    "first_doms.txt": "7d47c7c53e52b7bca0a1179b100bb8da",
    "results.csv": "b6d1aabdfc0f31dead6ac8541232a2c1",
}


def main(argv):
    p = argparse.ArgumentParser()
    p.add_argument("collection", help="catlas collection, e.g. twofoo")
    p.add_argument("filenames", nargs="+")
    p.add_argument("-q", "--quiet", action="store_true")
    args = p.parse_args()

    check_vals = saved_md5[args.collection]

    failed = False
    for filename in args.filenames:
        # do we know this file? if not, skip.
        basename = os.path.basename(filename)
        known_value = check_vals.get(basename)
        if not known_value:
            if not args.quiet:
                print("** not checking unknown file {}".format(basename))
            continue

        # calculate md5
        m = hashlib.md5()
        with open(filename, "rb") as fp:
            m.update(fp.read())

        digest = m.hexdigest()

        # check!
        if digest != known_value:
            failed = True
            print("previous {} != current {}".format(known_value, digest))

    if failed:
        return -1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
