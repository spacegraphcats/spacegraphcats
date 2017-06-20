#! /usr/bin/env python3

import argparse
import sys
from .walk_dbg import run
from .logging import log

# graph settings
DEFAULT_KSIZE = 31
DEFAULT_MEMORY = 1e8
DEFAULT_SCALED=100

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('seqfiles', nargs='+')
    p.add_argument('-o', '--output', default=None)
    p.add_argument('-k', '--ksize', default=DEFAULT_KSIZE, type=int)
    p.add_argument('-M', '--memory', default=DEFAULT_MEMORY,
                            type=float)
    p.add_argument('--force', action='store_true')
    p.add_argument('--label', action='store_true')
    p.add_argument('-l', '--loadgraph', type=str, default=None)
    p.add_argument('--no-assemble', help='do not output assembly',
                   action='store_true')
    args = p.parse_args()
    run(args)
    log(args.output, sys.argv)
