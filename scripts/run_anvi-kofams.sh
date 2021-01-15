#!/bin/bash -l

#SBATCH -J anvi_metabolism
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 28
#SBATCH --time=0-06:00:00
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=susheel.busi@uni.lu

echo "== Starting run at $(date)"

# Running ANVI'O kofams
work_dir="/work/projects/nomis/aaron_metaG/IMP3/saccharibacteria"
sample_file="saccharibacteria-genomes.txt"	# change as needed for genome of interest
threads=28

cd $work_dir
conda activate anvio-7

# run KOFAMS module
for i in `cat list`
do
    anvi-run-kegg-kofams -c "$i".db -T $threads
done 

# estimate metabolism for all genomes
anvi-estimate-metabolism -e $sample_file --matrix-format

echo "== Ending run at $(date)" 
