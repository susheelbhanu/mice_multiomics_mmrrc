##############################
# MODULES
import os, re
import glob
import pandas as pd
from functools import reduce

# ##############################
# # Parameters
# CORES=int(os.environ.get("CORES", 4))


##############################
# Paths
# SRC_DIR = srcdir("../scripts")
# ENV_DIR = srcdir("../envs")
# NOTES_DIR = srcdir("../notes")
# SUBMODULES= srcdir("../../submodules")
COV_DIR = config["cov_dir"]

##############################
# Dependencies 
PERL5LIB="/home/users/sbusi/perl5/lib/perl5${PERL5LIB:+:${PERL5LIB}}"
PERL_LOCAL_LIB_ROOT="/home/users/sbusi/perl5${PERL_LOCAL_LIB_ROOT:+:${PERL_LOCAL_LIB_ROOT}}"
PERL_MB_OPT="--install_base \"/home/users/sbusi/perl5\""
PERL_MM_OPT="INSTALL_BASE=/home/users/sbusi/perl5"


##############################
# default executable for snakemake
shell.executable("bash")


##############################
# working directory
workdir:
    config["work_dir"]


##############################
# Relevant directories
DATA_DIR = config["data_dir"]
DB_DIR = config["db_dir"]
FASTQ_DIR = config["fastq_dir"]
RESULTS_DIR = config["results_dir"]


##############################
# Steps
STEPS = config["steps"]


##############################
# Input
SAMPLES = [line.strip() for line in open("sample_list.txt").readlines()]
