import os
import yaml


def catlas_build(conf_file):
    "Produce the list of files output by 'spacegraphcats <config> build"
    with open(conf_file, "rt") as fp:
        jj = yaml.load(fp)

    catlas_base = jj["catlas_base"]
    ksize = jj["ksize"]
    radius = jj["radius"]

    dirname = "{}_k{}_r{}".format(catlas_base, ksize, radius)

    z = []
    z.append(os.path.join(catlas_base, "bcalm.{}.k31.unitigs.fa".format(catlas_base)))
    z.append(os.path.join(dirname, "catlas.csv"))
    z.append(os.path.join(dirname, "cdbg.gxt"))
    z.append(os.path.join(dirname, "contigs.fa.gz.indices"))
    z.append(os.path.join(dirname, "contigs.fa.gz.info.csv"))
    z.append(os.path.join(dirname, "contigs.fa.gz.mphf"))
    z.append(os.path.join(dirname, "first_doms.txt"))
    return z


def catlas_search(conf_file, cdbg_only=False, suffix=""):
    "Produce the list of files output by 'spacegraphcats <config> search"
    with open(conf_file, "rt") as fp:
        jj = yaml.load(fp)

    catlas_base = jj["catlas_base"]
    ksize = jj["ksize"]
    radius = jj["radius"]

    cdbg_str = ""
    if cdbg_only:
        cdbg_str = "_cdbg"
    dirname = "{}_k{}_r{}{}_search_oh0{}".format(
        catlas_base, ksize, radius, cdbg_str, suffix
    )

    filenames = jj["search"]
    z = []
    for x in filenames:
        x = os.path.basename(x)
        z.append(os.path.join(dirname, "{}.cdbg_ids.txt.gz".format(x)))
        z.append(os.path.join(dirname, "{}.contigs.sig".format(x)))
        z.append(os.path.join(dirname, "{}.frontier.txt.gz".format(x)))
        z.append(os.path.join(dirname, "{}.response.txt".format(x)))

    z.append(os.path.join(dirname, "results.csv"))

    return z


def catlas_extract(conf_file, cdbg_only=False, suffix=""):
    "Produce the list of files output by 'spacegraphcats <config> extract_contigs extract_reads"
    with open(conf_file, "rt") as fp:
        jj = yaml.load(fp)

    catlas_base = jj["catlas_base"]
    ksize = jj["ksize"]
    radius = jj["radius"]

    cdbg_str = ""
    if cdbg_only:
        cdbg_str = "_cdbg"
    dirname = "{}_k{}_r{}{}_search_oh0{}".format(
        catlas_base, ksize, radius, cdbg_str, suffix
    )

    filenames = jj["search"]
    z = []
    for x in filenames:
        x = os.path.basename(x)
        z.append(os.path.join(dirname, "{}.cdbg_ids.reads.gz".format(x)))
        z.append(os.path.join(dirname, "{}.cdbg_ids.contigs.fa.gz".format(x)))

    return z


def catlas_search_input(conf_file):
    "Produce the list of files output by 'spacegraphcats <config> search"
    with open(conf_file, "rt") as fp:
        jj = yaml.load(fp)

    filenames = jj["search"]
    return filenames
