import os
import yaml

def catlas_build(conf_file):
    "Produce the list of files output by 'conf/run <config> build"
    with open(conf_file, 'rt') as fp:
        jj = yaml.load(fp)

    catlas_base = jj['catlas_base']
    ksize = jj['ksize']
    radius = jj['radius']

    dirname = '{}_k{}_r{}'.format(catlas_base, ksize, radius)

    z = []
    z.append(catlas_base + '/bcalm.{}.k31.unitigs.fa'.format(catlas_base))
    z.append(dirname + '/catlas.csv')
    z.append(dirname + '/cdbg.gxt')
    z.append(dirname + '/contigs.fa.gz')
    z.append(dirname + '/contigs.fa.gz.indices')
    z.append(dirname + '/contigs.fa.gz.info.csv')
    z.append(dirname + '/contigs.fa.gz.mphf')
    z.append(dirname + '/first_doms.txt')
    return z


def catlas_search(conf_file, cdbg_only=False, suffix=''):
    "Produce the list of files output by 'conf/run <config> search"
    with open(conf_file, 'rt') as fp:
        jj = yaml.load(fp)

    catlas_base = jj['catlas_base']
    ksize = jj['ksize']
    radius = jj['radius']

    cdbg_str = ""
    if cdbg_only:
        cdbg_str = '_cdbg'
    dirname = '{}_k{}_r{}{}_search_oh0{}'.format(catlas_base, ksize, radius,
                                                 cdbg_str, suffix)

    filenames = jj['search']
    z = []
    for x in filenames:
        x = os.path.basename(x)
        z.append(dirname + '/{}.cdbg_ids.txt.gz'.format(x))
        z.append(dirname + '/{}.contigs.sig'.format(x))
        z.append(dirname + '/{}.frontier.txt.gz'.format(x))
        z.append(dirname + '/{}.response.txt'.format(x))

    z.append(dirname + '/results.csv')

    return z


def catlas_extract(conf_file, cdbg_only=False, suffix=''):
    "Produce the list of files output by 'conf/run <config> extract_contigs extract_reads"
    with open(conf_file, 'rt') as fp:
        jj = yaml.load(fp)

    catlas_base = jj['catlas_base']
    ksize = jj['ksize']
    radius = jj['radius']

    cdbg_str = ""
    if cdbg_only:
        cdbg_str = '_cdbg'
    dirname = '{}_k{}_r{}{}_search_oh0{}'.format(catlas_base, ksize, radius,
                                                 cdbg_str, suffix)

    filenames = jj['search']
    z = []
    for x in filenames:
        x = os.path.basename(x)
        z.append(dirname + '/{}.cdbg_ids.reads.fa.gz'.format(x))
        z.append(dirname + '/{}.cdbg_ids.contigs.fa.gz'.format(x))

    return z


def catlas_search_input(conf_file):
    "Produce the list of files output by 'conf/run <config> search"
    with open(conf_file, 'rt') as fp:
        jj = yaml.load(fp)

    filenames = jj['search']
    return filenames
