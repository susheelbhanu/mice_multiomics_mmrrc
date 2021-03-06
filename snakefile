"""
Author: Susheel Bhanu BUSI
Affiliation: ESB group LCSB UniLU
Date: [2019-06-02]
Run: snakemake -s snakefile
Latest modification:
"""

import os
import glob
# import pandas 

configfile:"config.yaml"
DATA_DIR=config['data_dir']
RESULTS_DIR=config['results_dir']
SAMPLES=[line.strip() for line in open("sample_list", 'r')]	# if using a sample list instead of putting them in a config file
# SAMPLES=config['samples']

# print(f"number of cores: {workflow.cores}")

###########
rule all: 
		input: 
			expand(os.path.join(DATA_DIR, "IMP3/{sample}/{sample}_{type}_R{read}.fastq.gz"), type=["DNA", "RNA"], read=["1", "2"], sample=SAMPLES),
			expand(os.path.join(DATA_DIR, "IMP3/{sample}/{sample}_config.yaml"), sample=SAMPLES),
			expand(os.path.join(DATA_DIR, "IMP3/{sample}/launchIMP.sh"), sample=SAMPLES),
			expand(os.path.join(DATA_DIR, "IMP3/{sample}/{sample}.runIMP.sh"), sample=SAMPLES),
#			expand(os.path.join(DATA_DIR, "IMP3/{sample}/run1/status/all.done"), sample=SAMPLES)
#			expand(os.path.join(DATA_DIR, "IMP3/{sample}/run1/Analysis/eukrep/{sample}_eukaryotic_contigs.fa"), sample=SAMPLES),
#			expand(os.path.join(DATA_DIR, "IMP3/{sample}/run1/Analysis/nonpareil/{sample}.fasta"), sample=SAMPLES),
#			expand(os.path.join(DATA_DIR, "IMP3/{sample}/run1/Analysis/metaxa2/{sample}_output/{sample}_metaxa_output.taxonomy.txt"), sample=SAMPLES),
#			expand(os.path.join(DATA_DIR, "IMP3/{sample}/run1/Analysis/metaxa2/{sample}_output/LEVEL-7_{sample}_rarefaction_out"), sample=SAMPLES)
#			expand(os.path.join(DATA_DIR, "IMP3/{sample}/run1/Assembly/intermediary/mg.megahit_preprocessed.1/options.json"), sample=SAMPLES)

################################
# rules for files and analyses #
################################
rule folders:
		input:
			in1=os.path.join(DATA_DIR, "IMP3/{sample}_DNA_R1.fastq.gz"),
			in2=os.path.join(DATA_DIR, "IMP3/{sample}_DNA_R2.fastq.gz"),
			in3=os.path.join(DATA_DIR, "IMP3/{sample}_RNA_R1.fastq.gz"),
			in4=os.path.join(DATA_DIR, "IMP3/{sample}_RNA_R2.fastq.gz")
		output:	
			fout1=os.path.join(DATA_DIR, "IMP3/{sample}/{sample}_DNA_R1.fastq.gz"),
			fout2=os.path.join(DATA_DIR, "IMP3/{sample}/{sample}_DNA_R2.fastq.gz"),
			fout3=os.path.join(DATA_DIR, "IMP3/{sample}/{sample}_RNA_R1.fastq.gz"),
			fout4=os.path.join(DATA_DIR, "IMP3/{sample}/{sample}_RNA_R2.fastq.gz"),
		shell:	"""
			date &&\
			mv -v {input.in1} {output.fout1} && mv -v {input.in2} {output.fout2} &&\
			mv -v {input.in3} {output.fout3} && mv -v {input.in4} {output.fout4} &&\
			date
			"""

rule imp_files:
		output:
			tout1=os.path.join(DATA_DIR, "IMP3/{sample}/{sample}_config.yaml"),
                        tout2=os.path.join(DATA_DIR, "IMP3/{sample}/launchIMP.sh"),
                        tout3=os.path.join(DATA_DIR, "IMP3/{sample}/{sample}.runIMP.sh")
		shell:	"""
			date &&\
			cp -v IMP_config.yaml {output.tout1} && cp -v launchIMP.sh {output.tout2} && cp -v runIMP.sh {output.tout3} &&\
			sed -i 's/\"\$sample\"/{wildcards.sample}/g' {output.tout1} &&\
			sed -i 's/\"\$sample\"/{wildcards.sample}/g' {output.tout2} &&\
			sed -i 's/\"\$sample\"/{wildcards.sample}/g' {output.tout3} &&\
			date
			"""

#rule imp_launch:
#		input:
#			l1=os.path.join(DATA_DIR, "IMP3/{sample}")
#		output:
#			lout1=os.path.join(DATA_DIR, "IMP3/{sample}/run1/status/all.done"),
#			lout2=os.path.join(DATA_DIR, "IMP3/{sample}/run1/Assembly/mg.assembly.merged.fa"),
#			lout3=os.path.join(DATA_DIR, "IMP3/{sample}/run1/Preprocessing/mg.r1.trimmed.hg38_filtered.fq"),
#			lout4=os.path.join(DATA_DIR, "IMP3/{sample}/run1/Preprocessing/mg.r2.trimmed.hg38_filtered.fq")
#		threads:config["threads"]
#		shell:	"""
#			date
#			chmod -R 775 {input.l1} && cd {input.l1} 
#			IMPROOT="/work/projects/ecosystem_biology/local_tools/IMP3"
#			snakemake -s $IMPROOT/Snakefile --configfile {wildcards.sample}_config.yaml --cores {threads} --use-conda --conda-prefix $IMPROOT/conda --unlock
#			snakemake -s $IMPROOT/Snakefile --configfile {wildcards.sample}_config.yaml --cores {threads} --use-conda --conda-prefix $IMPROOT/conda
#			date
#			"""

#rule decompress:
#		input:	os.path.join(DATA_DIR, "IMP3/{sample}/run1/Assembly/intermediary.tar.gz")
#		output:	os.path.join(DATA_DIR, "IMP3/{sample}/run1/Assembly/intermediary/mg.megahit_preprocessed.1/options.json")
#		shell:	"""
#			date &&\
#			tar -xzvf {input} -C $(dirname {input}) &&\
#			date
#			"""
#
#rule copy_updated_config:
#		output:	os.path.join(DATA_DIR, "IMP3/{sample}/{sample}_config.yaml")
#		shell:	"""
#			date &&\
#			cp -v intermediary_IMP_config.yaml {output} &&\
#			sed -i 's/\"\$sample\"/{wildcards.sample}/g' {output} &&\
#			date
#			"""
			

############
# Analyses #
############
rule Eukrep_fasta:
		input:
			pru1=os.path.join(DATA_DIR, "IMP3/{sample}/run1/Assembly/mg.assembly.merged.fa")
		threads:config["threads"]
		output:
			prout=os.path.join(DATA_DIR, "IMP3/{sample}/run1/Analysis/eukrep/{sample}_eukaryotic_contigs.fa")
		shell:	"""
			source /home/users/sbusi/apps/Eukrep/bin/activate
			date
			EukRep -i {input.pru1} -o {output.prout} --min 2000 -m strict
			date
			"""

#####################################
##### NONPAREIL & METAXA2 ###########
#####################################                        
rule nonpareil:
		input:
			np1=os.path.join(DATA_DIR, "IMP3/{sample}/run1/Preprocessing/mg.r1.trimmed.hg38_filtered.fq")
		threads:config["threads"]
		conda:	"/home/users/sbusi/apps/environments/conda-env.yml"
		output:
			np2=os.path.join(DATA_DIR, "IMP3/{sample}/run1/Analysis/nonpareil/{sample}.fasta")
		shell:	"""
			date
			cat {input.np1} | paste - - - - | awk 'BEGIN{{{{FS="\t"}}{{}}{{print ">"substr($1,2)"\n"$2}}}}' > {output.np2}
			nonpareil -s {output.np2} -T kmer -f fasta -b {wildcards.sample}
			date
			"""

rule metaxa2:
		input:
			mt1=os.path.join(DATA_DIR, "IMP3/{sample}/run1/Preprocessing/mg.r1.trimmed.hg38_filtered.fq"),
			mt2=os.path.join(DATA_DIR, "IMP3/{sample}/run1/Preprocessing/mg.r2.trimmed.hg38_filtered.fq")
		threads:config["threads"]
		conda:	"/home/users/sbusi/apps/environments/conda_env.yml"
		output:
#			metout=directory(os.path.join(DATA_DIR, "IMP3/{sample}/run1/Analysis/metaxa2/{sample}_output")),
			metout=os.path.join(DATA_DIR, "IMP3/{sample}/run1/Analysis/metaxa2/{sample}_output/{sample}_metaxa_output.taxonomy.txt")
		shell:	"""
			date
			metaxa2 --threads {threads} -f fastq -1 {input.mt1} -2 {input.mt2} -o $(dirname {output.metout})
			date
			"""

rule rarefaction_metaxa2:
		input:
			ra1=os.path.join(DATA_DIR, "IMP3/{sample}/run1/Analysis/metaxa2/{sample}_output/{sample}_metaxa_output.taxonomy.txt")
		threads:config["threads"]
		conda:	"/home/users/sbusi/apps/environments/conda_env.yml"
		output:
			ra2=os.path.join(DATA_DIR, "IMP3/{sample}/run1/Analysis/metaxa2/{sample}_output/LEVEL-7_{sample}_rarefaction_out")
		shell:	"""
			date
			metaxa2_rf --threads {threads} -n 7 --resamples 10000 --scale 1000000 -i {input.ra1} -o {output.ra2}
			date
			"""


