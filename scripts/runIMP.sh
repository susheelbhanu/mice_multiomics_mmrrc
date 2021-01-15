#!/bin/bash -l

export PATH=$PATH:/work/projects/ecosystem_biology/local_tools/bin
export PATH=/work/projects/ecosystem_biology/local_tools/minicondaESB/condabin:/work/projects/ecosystem_biology/local_tools/bin:$PATH

IMPROOT=/work/projects/ecosystem_biology/local_tools/IMP3
THREADS=8

snakemake -s $IMPROOT/Snakefile --configfile "$sample"_config.yaml -j $THREADS --use-conda --conda-prefix $IMPROOT/conda --unlock
snakemake -s $IMPROOT/Snakefile --configfile "$sample"_config.yaml -j $THREADS --use-conda --conda-prefix $IMPROOT/conda
# snakemake -s $IMPROOT/Snakefile --configfile "$sample".config.json -j $THREADS --use-conda --conda-prefix $IMPROOT/conda --report "$sample".report.html

# moving the files to a different location
# rsync --remove-source-files -auv --no-p --no-g -P /scratch/users/sbusi/aaron_metaG/IMP3/"$sample" /work/projects/nomis/aaron_metaG/IMP3/.

