# eukaryotes workflow
# Runs EUKulele to obtain eukaryotic abundances

import os
import pandas as pd

# CONFIG
# config file
configfile: "config.yaml"

# default executable for snakmake
shell.executable("bash")

# Paths
SRC_DIR = srcdir("scripts")
ENV_DIR = srcdir("envs")

# Parameters
DATA_DIR=config["data_dir"]
RESULTS_DIR=config["results_dir"]
# SAMPLES=[line.strip() for line in open("list", 'r')]
SAMPLES=config["samples"]

localrules: merge_cov, merge_cov_all, euk_abundances, euk_abundances_all

###########################
# default
rule eukaryotes: 
    input: 
        os.path.join(RESULTS_DIR, "eukulele_output")
#        expand(os.path.join(RESULTS_DIR, "euk_tax_cov/{sample}_eukaryotes.txt"), sample=SAMPLES),
#        os.path.join(RESULTS_DIR, "euk_tax_cov/merged_eukaryote_abundances.txt"),
#        os.path.join(RESULTS_DIR, "euk_tax_cov/merged_eukulele_all_abundances.txt")
    output:
        touch("status/eukaryotes.done")


###################################
# rules for Eukaryotic annotation #
###################################
rule file_prep:
    input:
        os.path.join(DATA_DIR, "{sample}/run1/Analysis/annotation/prokka.faa")
    output:
        os.path.join(RESULTS_DIR, "eukulele_input/{sample}.faa") 
    log:
        os.path.join(RESULTS_DIR, "logs/eukulele_{sample}_prep.log")
    shell:
        "(date && ln -vs {input} {output} && date) &> >(tee {log})"

rule eukulele:
    input:    
        expand(os.path.join(RESULTS_DIR, "eukulele_input/{sample}.faa"), sample=[open('sample_list').readline().strip("\n")])
    output:    
        directory(os.path.join(RESULTS_DIR, "eukulele_output"))
    threads:
        config['threads']
    conda:
        os.path.join(ENV_DIR, "eukulele.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/eukulele.log")
    message:
        "Running Eukaryotic annotation with Eukulele"
    shell:
        "(date && EUKulele --sample_dir $(dirname {input}) -o {output[0]} -m mets && date) &> >(tee {log})"

rule merge_cov:
    input:
        cov=os.path.join(DATA_DIR, "{sample}/run1/Stats/mg/annotation/mg.gene_depth.avg")
    output:
        os.path.join(RESULTS_DIR, "euk_tax_cov/{sample}_eukaryotes.txt")
    params:
        tax=os.path.join(RESULTS_DIR, "eukulele_output/taxonomy_estimation/{sample}-estimated-taxonomy.out")
    log:
        os.path.join(RESULTS_DIR, "logs/cov_merge.{sample}.log")
    message:
        "Merging Eukaryotic taxonomy with coverage for {wildcards.sample}"
    run:
        tax=pd.read_csv(params.tax, header=0, sep="\t", index_col=0)
        cov=pd.read_csv(input.cov, header=None, sep="\t")
        cov.columns=['transcript_name', 'coverage']
        
        # keeping only rows that contain 'Eukary' in the 'full_classification' column
        euks=tax[tax['full_classification'].str.contains("Eukary", na=False)].drop(['counts'], axis=1)
        # keeping only rows that have at least 70% 'max_pid'
        filt_euks=euks.query("max_pid >=70")

        # merging taxonomy with coverage
        merged=filt_euks.merge(cov, how="left", on="transcript_name")
        filt_merged=merged[['full_classification','classification', 'coverage']]

        # Grouping same taxonomy and getting sum of coverage
        final=filt_merged.groupby(['full_classification','classification'], as_index=False)['coverage'].sum()
        
        # writing to file
        final.to_csv(output[0], index=None, sep="\t")

rule euk_abundances:
    input:
        expand(os.path.join(RESULTS_DIR, "euk_tax_cov/{sample}_eukaryotes.txt"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "euk_tax_cov/merged_eukaryote_abundances.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/euk_abundances.log")
    message:
        "Merging all eukaryote abundances for ROCK samples"
    run:
        # Collecting all files in folder 
        directory=os.path.dirname(input[0])
        os.chdir(directory)
	
        # verify the path using getcwd() 
        cwd = os.getcwd() 
  
        # print the current directory 
        print("Current working directory is:", cwd) 

        mylist=[f for f in glob.glob("*.txt")]
        mylist

        # making individual dataframes for each file
        dataframes= [ pd.read_csv( f, header=0, sep="\t", usecols=["full_classification", "coverage"]) for f in mylist ] # add arguments as necessary to the read_csv method

        # Merging all files based on common column
        merged=reduce(lambda left,right: pd.merge(left,right,on='full_classification', how='outer'), dataframes)

        # Giving appropriate column names
        names=['full_classification']+mylist
        new_cols=list(map(lambda x: x.replace('_eukaryotes.txt',''),names))
        merged.columns=new_cols

        # checking if any values are "NA"
        merged.isnull().values.any()
        # if "NA" run the following
        merged.fillna('', inplace=True)

        # Removing rows with all zeroes (0 or 0.0)
        merged.set_index('full_classification', inplace=True)  # first to make first column as rownames
        edited=merged.loc[~(merged==0).all(axis=1)]

        # Writing file without zeroes
        edited.to_csv(output[0], sep='\t', index=True, header=True)

rule merge_cov_all:
    input:
        cov=os.path.join(DATA_DIR, "{sample}/run1/Stats/mg/annotation/mg.gene_depth.avg")
    output:
        os.path.join(RESULTS_DIR, "euk_tax_cov/{sample}_eukulele_ALL.txt")
    params:
        tax=os.path.join(RESULTS_DIR, "eukulele_output/taxonomy_estimation/{sample}-estimated-taxonomy.out")
    log:
        os.path.join(RESULTS_DIR, "logs/cov_merge_all.{sample}.log")
    message:
        "Merging Eukaryotic taxonomy with coverage for {wildcards.sample}"
    run:
        tax=pd.read_csv(input.tax, header=0, sep="\t", index_col=0)
        cov=pd.read_csv(input.cov, header=None, sep="\t")
        cov.columns=['transcript_name', 'coverage']

        # keeping only rows that contain 'Eukary' in the 'full_classification' column
        euks=tax.drop(['counts'], axis=1)

        # merging taxonomy with coverage
        merged=euks.merge(cov, how="left", on="transcript_name")
        filt_merged=merged[['full_classification','classification', 'coverage']]

        # Grouping same taxonomy and getting sum of coverage
        final=filt_merged.groupby(['full_classification','classification'], as_index=False)['coverage'].sum()

        # writing to file
        final.to_csv(output[0], index=None, sep="\t")

rule euk_abundances_all:
    input:
        expand(os.path.join(RESULTS_DIR, "euk_tax_cov/{sample}_eukulele_ALL.txt"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "euk_tax_cov/merged_eukulele_all_abundances.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/euk_all_abundances.log")
    message:
        "Merging all eukaryote abundances for ROCK samples"
    run:
        # Collecting all files in folder
        directory=os.path.dirname(input[0])
        os.chdir(directory)

        # verify the path using getcwd()
        cwd = os.getcwd()

        # print the current directory
        print("Current working directory is:", cwd)

        mylist=[f for f in glob.glob("*ALL.txt")]
        mylist

        # making individual dataframes for each file
        dataframes= [ pd.read_csv( f, sep="\t", usecols=['full_classification', 'coverage']) for f in mylist ] # add arguments as necessary to the read_csv method

        # Merging all files based on common column
        merged=reduce(lambda left,right: pd.merge(left,right,on='full_classification', how='outer'), dataframes)

        # Giving appropriate column names
        names=['full_classification']+mylist
        new_cols=list(map(lambda x: x.replace('_eukulele_all.txt',''),names))
        merged.columns=new_cols

        # checking if any values are "NA"
        merged.isnull().values.any()
        # if "NA" run the following
        merged.fillna('', inplace=True)

        # Removing rows with all zeroes (0 or 0.0)
        merged.set_index('full_classification', inplace=True)  # first to make first column as rownames
        edited=merged.loc[~(merged==0).all(axis=1)]

        # Writing file without zeroes
        edited.to_csv(output[0], sep='\t', index=True, header=True)
