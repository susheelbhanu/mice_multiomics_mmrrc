#!/bin/bash -l
sbatch --time=120:00:0 -N1 -n8 -p bigmem -q long  -J "$sample" ./"$sample".runIMP.sh
