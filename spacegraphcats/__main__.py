"""
Execution script for snakemake workflows for spacegraphcats.
"""
import sys
import argparse
import os.path
import snakemake
import sys
import pprint


thisdir = os.path.abspath(os.path.dirname(__file__))


def main():
    parser = argparse.ArgumentParser(description='run snakemake workflows for spacegraphcats', usage='''run <configfile.yaml> [<target1> ...]

Run snakemake workflows for spacegraphcats, using the given config file.

Targets:

   build            - builds catlas (default)
   searchquick      - do a quick search (only a few files)
   search           - do a full search for this data set (many files)
   extract_reads    - extract reads for search results (many files)
   extract_contigs  - extract contigs for search results (many files)
   clean            - remove the primary catlas build files

For a quickstart, run this:

   conf/run dory-test searchquick

from the main spacegraphcats directory.
.
''')

    parser.add_argument('configfile')
    parser.add_argument('targets', nargs='*', default=['build'])
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('--nolock', action='store_true')
    parser.add_argument('--overhead', type=float, default=None)
    parser.add_argument('--experiment', default=None)
    parser.add_argument('--radius', type=int, default=None)
    parser.add_argument('--cdbg-only', action='store_true',
                        help='for paper evaluation purposes')
    args = parser.parse_args()

    # first, find the Snakefile
    snakefile = os.path.join(thisdir, './conf', 'Snakefile')
    if not os.path.exists(snakefile):
        sys.stderr.write('Error: cannot find Snakefile at {}\n'.format(snakefile))
        sys.exit(-1)

    # next, find the config file
    configfile = None
    if os.path.exists(args.configfile) and not os.path.isdir(args.configfile):
        configfile = args.configfile
    else:
        for suffix in ('', '.json', '.yaml'):
            tryfile = os.path.join(thisdir, './conf', args.configfile + suffix)
            if os.path.exists(tryfile) and not os.path.isdir(tryfile):
                sys.stderr.write('Found configfile at {}\n'.format(tryfile))
                configfile = tryfile
                break

    if not configfile:
        sys.stderr.write('Error: cannot find configfile {}\n'.format(args.configfile))
        sys.exit(-1)

    # build config override dict
    config = dict()
    if args.overhead is not None:
        config['overhead'] = args.overhead
    if args.experiment is not None:
        config['experiment'] = args.experiment
    if args.radius is not None:
        config['radius'] = args.radius
    if args.cdbg_only:
        config['cdbg_only'] = True

    print('--------')
    print('details!')
    print('\tsnakefile: {}'.format(snakefile))
    print('\tconfig: {}'.format(configfile))
    print('\ttargets: {}'.format(repr(args.targets)))
    if config:
        print('\toverride: {}'.format(pprint.pformat(config)))
    print('--------')

    # run!!
    status = snakemake.snakemake(snakefile, configfile=configfile,
                                 targets=args.targets, printshellcmds=True,
                                 dryrun=args.dry_run,
                                 lock=not args.nolock, config=config,
                                 verbose=args.verbose, debug_dag=args.debug)

    if status: # translate "success" into shell exit code of 0
       return 0
    return 1


if __name__ == '__main__':
    sys.exit(main())
