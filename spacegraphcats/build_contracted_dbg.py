import argparse
from walk_dbg import run

# graph settings
DEFAULT_KSIZE = 31
DEFAULT_MEMORY = 1e8

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
    args = p.parse_args()

    run(args)
