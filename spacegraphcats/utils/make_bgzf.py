#! /usr/bin/env python
"""
Convert an input file into a BGZF file that supports random indexing by offset.
"""
import screed
from .bgzf import bgzf
import argparse
import sys


def main(argv=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument("input_files", nargs="+")
    p.add_argument("-o", "--output")
    args = p.parse_args(argv)

    assert args.output, "must specify -o"

    output_filename = args.output
    outfp = bgzf.open(output_filename, "wb")

    print("output file will be {}".format(output_filename))

    for input_file in args.input_files:
        print("turning {} into a block-gzipped (BGZF) file".format(input_file))
        for n, record in enumerate(screed.open(input_file)):
            offset = outfp.tell()
            if hasattr(record, "quality"):
                outfp.write(
                    "@{}\n{}\n+\n{}\n".format(
                        record.name, record.sequence, record.quality
                    )
                )
            else:
                outfp.write(">{}\n{}\n".format(record.name, record.sequence))
            if n % 100000 == 0:
                print("offset for {} is {}".format(n, offset), end="\r")
        print("")

    outfp.close()

    print("done!")

    return 0


if __name__ == "__main__":
    sys.exit(main())
