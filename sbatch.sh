#!/bin/bash -l

##############################
# SLURM
# NOTE: used for this script only, NOT for the snakemake call below

#SBATCH -J AE_euk
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=4-00:00:00
#SBATCH -p batch
#SBATCH --qos=qos-batch

##############################
# SNAKEMAKE

# conda env name
SMK_ENV="snakemake" # USER INPUT REQUIRED
# number of cores for snakemake
SMK_CORES=6
# number of jobs for snakemake
SMK_JOBS=6
# snakemake file
SMK_SMK="eukaryote_snakefile.smk"
# config file
SMK_CONFIG="config/config.yaml" # USER INPUT REQUIRED
# slurm config file
SMK_SLURM="config/slurm.yaml"
# slurm cluster call
SMK_CLUSTER="sbatch -p {cluster.partition} -q {cluster.qos} {cluster.explicit} -N {cluster.nodes} -n {cluster.n} -c {threads} -t {cluster.time} --job-name={cluster.job-name}"


##############################
# IMP

# activate the env
conda activate ${SMK_ENV}

# run the pipeline
snakemake -s ${SMK_SMK} --cores ${SMK_CORES} --jobs ${SMK_JOBS} --local-cores 1 \
--configfile ${SMK_CONFIG} --use-conda --conda-prefix ${CONDA_PREFIX}/pipeline \
--cluster-config ${SMK_SLURM} --cluster "${SMK_CLUSTER}" --rerun-incomplete -rp

