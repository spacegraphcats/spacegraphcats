# Integrating spacegraphcats into a snakemake workflow 

This page is under construction, please check back as we add more documentation!
Snakemake is a workflow manager that automates and scales bioinformatics pipelines, and makes them portable across computing environments.

## A sample snakemake workflow

## Configuration files

### Generating configuration files prior to running snakemake

### Using snakemake to automatically generate configuration files

In some workflows, the queries may not be pre-defined but instead may be first produced by the workflow.
In this case, [snakemake checkpoints](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution) can be used to automatically and flexibly generate spacegraphcats configuration files based on the results of the workflow.
