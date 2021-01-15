#!/bin/bash -l

export PYTHONNOUSERSITE=TRUE

SLURM_ARGS="-p {cluster.partition} --qos {cluster.qos} -N 1 -n {cluster.n} -c {cluster.ncpus} -t {cluster.time} --job-name={cluster.job-name} -o {cluster.output} -e {cluster.error}"
(date; conda activate snakemake; snakemake -p -s snakefile --unlock; snakemake -j 8 --latency-wait 30 --cluster-config cluster.json --cluster "sbatch $SLURM_ARGS" -s snakefile -p --use-conda ; date) &> smk.log
