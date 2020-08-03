"Enable python -m spacegraphcats.click"
import sys
import os
import subprocess
import glob

import click

def get_snakefile_path(name):
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, 'conf', name)
    return snakefile


def get_package_configfile(filename):
    thisdir = os.path.dirname(__file__)
    configfile = os.path.join(thisdir, 'conf', filename)
    return configfile


def run_snakemake(configfile, no_use_conda=False, verbose=False,
                  snakefile_name='Snakefile', outdir=None, extra_args=[]):
    # find the Snakefile relative to package path
    snakefile = get_snakefile_path(snakefile_name)

    # basic command
    cmd = ["snakemake", "-s", snakefile]

    # add --use-conda
    if not no_use_conda:
        cmd += ["--use-conda"]

    # add --outdir
    if outdir:
        cmd += ["--config", f"outdir={outdir}"]

    # snakemake sometimes seems to want a default -j; set it to 1 for now.
    # can overridden later on command line.
    cmd += ["-j", "1"]

    # add rest of snakemake arguments
    cmd += list(extra_args)

    # add configfile - try looking for it a few different ways.
    configfiles = []
    if os.path.isfile(configfile):
        configfiles = [configfile]
    elif os.path.isfile(get_package_configfile(configfile)):
        configfiles = [get_package_configfile(configfile)]
    else:
        for suffix in '.yaml', '.conf':
            tryfile = configfile + suffix
            if os.path.isfile(tryfile):
                configfiles = [tryfile]
                break

            tryfile = get_package_configfile(tryfile)
            if os.path.isfile(tryfile):
                configfiles = [tryfile]
                break

    if not configfiles:
        raise ValueError(f"cannot find config file '{configfile}'")

    cmd += ["--configfile"] + configfiles

    if verbose:
        print('final command:', cmd)

    # runme
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f'Error in snakemake invocation: {e}', file=sys.stderr)
        return e.returncode

    return 0

#
# actual command line functions
#

@click.group()
def cli():
    pass

# create a run subcommand that by default passes all of its arguments
# on to snakemake (after setting Snakefile and config)
@click.command(context_settings={"ignore_unknown_options": True})
@click.argument('configfile')
@click.option('--no-use-conda', is_flag=True, default=False)
@click.option('--verbose', is_flag=True)
@click.option('--outdir', nargs=1)
@click.argument('snakemake_args', nargs=-1)
def run(configfile, snakemake_args, no_use_conda, verbose, outdir):
    "execute spacegraphcats workflow (using snakemake underneath)"
    run_snakemake(configfile, snakefile_name='Snakefile',
                  no_use_conda=no_use_conda, verbose=verbose,
                  extra_args=snakemake_args, outdir=outdir)

# build command -- run the snakemake target build.
@click.command(context_settings={"ignore_unknown_options": True})
@click.argument('configfile')
@click.option('--no-use-conda', is_flag=True, default=False)
@click.option('--verbose', is_flag=True)
@click.option('--outdir', nargs=1)
@click.argument('snakemake_args', nargs=-1)
def build(configfile, snakemake_args, no_use_conda, verbose, outdir):
    "build spacegraphcats indices"
    run_snakemake(configfile, snakefile_name='Snakefile',
                  no_use_conda=no_use_conda, verbose=verbose,
                  extra_args=['build'] + list(snakemake_args), outdir=outdir)

# search command -- run the snakemake target search.
@click.command(context_settings={"ignore_unknown_options": True})
@click.argument('configfile')
@click.option('--no-use-conda', is_flag=True, default=False)
@click.option('--verbose', is_flag=True)
@click.option('--outdir', nargs=1)
@click.argument('snakemake_args', nargs=-1)
def search(configfile, snakemake_args, no_use_conda, verbose, outdir):
    "search spacegraphcats indices with a query"
    run_snakemake(configfile, snakefile_name='Snakefile',
                  no_use_conda=no_use_conda, verbose=verbose,
                  extra_args=['search'] + list(snakemake_args), outdir=outdir)

# 'check' command
@click.command()
@click.argument('configfile')
def check(configfile):
    "check configuration"
    run_snakemake(configfile, extra_args=['check'])

# 'showconf' command
@click.command()
@click.argument('configfile')
def showconf(configfile):
    "show full configuration"
    run_snakemake(configfile, extra_args=['showconf'])

# 'info' command
@click.command()
def info():
    "provide basic install/config file info"
    from .version import version
    print(f"""
This is spacegraphcats version v{version}

Package install path: {os.path.dirname(__file__)}
snakemake Snakefile: {get_snakefile_path('Snakefile')}
""")

# 'init' command
@click.command()
@click.argument('configfile')
@click.option('-f', '--force', is_flag=True)
def init(configfile, force):
    "create a new, empty config file."
    stubname = os.path.basename(configfile)
    if configfile.endswith('.conf'):
        stubname = stubname[:-5]
    else:
        configfile += '.conf'

    if os.path.exists(configfile) and not force:
        print(f"** ERROR: configfile '{configfile}' already exists.")
        return -1

    # @CTB
    print(f"creating configfile '{configfile}' for project '{stubname}'")
    with open(configfile, 'wt') as fp:
        fp.write(\
f"""\
# basic configuration:
catlas_base: twofoo
input_sequences:
- twofoo.fq.gz
ksize: 31
radius: 1
search:
- data/2.fa.gz
- data/47.fa.gz
- data/63.fa.gz
searchquick:
- data/2.fa.gz

## advanced features
# multifasta query
multifasta_query_sig: data/63-os223.sig
multifasta_reference:
- data/twofoo-genes.fa.gz
multifasta_scaled: 1000

# hashval queries
hashval_ksize: 51
hashval_queries: data/twofoo-k51-hashval-queries.txt

# shadow ratio foo
shadow_ratio_maxsize: 1000
# END
""")

cli.add_command(run)
cli.add_command(build)
cli.add_command(search)
cli.add_command(check)
cli.add_command(showconf)
cli.add_command(info)
cli.add_command(init)

def main():
    cli()

if __name__ == '__main__':
    main()
