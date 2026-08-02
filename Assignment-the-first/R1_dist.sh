#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --job-name=R1_dist

DATA=/projects/bgmp/shared/2017_sequencing/
R1=$DATA/1294_S1_L008_R1_001.fastq.gz
R2=$DATA/1294_S1_L008_R2_001.fastq.gz
R3=$DATA/1294_S1_L008_R3_001.fastq.gz
R4=$DATA/1294_S1_L008_R4_001.fastq.gz

/usr/bin/time -v python bioreads_dist.py -f $R1 -r 1
