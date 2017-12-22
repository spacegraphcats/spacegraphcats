# pull values in from config file:
catlas_dir=config['catlas_dir']
DATASETS=config['DATASETS']
ksize=config['ksize']
radius=config['radius']
seeds=config['seeds']

# internal definitions for convenience:
python=sys.executable  # define as the version of Python running snakemake

# rewrite the datasets into bcalm format
def write_bcalm_in(datasets):
    x = []
    for d in datasets:
        x.append('-in {}'.format(d))
    return " ".join(x)

BCALM_INPUT=write_bcalm_in(DATASETS)

### rules!

# build catlas & minhashes needed for search
rule all:
    input:
        expand("{catlas_dir}/catlas.csv", catlas_dir=catlas_dir),
        expand("{catlas_dir}/minhashes.db.k{ksize}.s1000.abund0.seed{seed}", catlas_dir=catlas_dir, seed=seeds, ksize=ksize)

# build cDBG using bcalm
rule bcalm_cdbg:
     input:
        DATASETS
     output:
        "bcalm.{catlas_dir}.k{ksize}.unitigs.fa"
     shell:
        # @CTB here we run into the problem that bcalm wants
        # '-in file1 -in file2', so I am using 'params' and 
        "bcalm {BCALM_INPUT} -out bcalm.{wildcards.catlas_dir}.k{wildcards.ksize} -kmer-size {wildcards.ksize} >& {output}.log.txt"

# build catlas input from bcalm output by reformatting
rule bcalm_catlas_input:
     input:
        expand("bcalm.{{catlas_dir}}.k{ksize}.unitigs.fa", ksize=ksize)
     output:
        "{catlas_dir}/cdbg.gxt",
        "{catlas_dir}/contigs.fa.gz"
     shell:
        "{python} bcalm-to-gxt.py {input} {output}"

# build catlas!
rule build_catlas:
     input:
        "{catlas_dir}/cdbg.gxt",
        "{catlas_dir}/contigs.fa.gz",
     output:
        "{catlas_dir}/first_doms.txt",
        "{catlas_dir}/catlas.csv"
     shell:
        "{python} -m spacegraphcats.catlas {catlas_dir} {radius}"

# build minhash databases
rule minhash_db:
     input:
        "{catlas_dir}/cdbg.gxt",
        "{catlas_dir}/contigs.fa.gz",
        "{catlas_dir}/first_doms.txt",
     output:
        "{catlas_dir}/minhashes.db.k{ksize}.s1000.abund0.seed{seed}"
     shell:
        "{python} -m search.make_catlas_minhashes -k {wildcards.ksize} --seed={wildcards.seed} --scaled=1000 {catlas_dir}"
