#!/bin/bash

#SBATCH -J enrich 
#SBATCH -n 1                # Use n cores for one job 
#SBATCH -t 0-12:00                # Runtime in D-HH:MM 
#SBATCH -p serial_requeue                # Partition to submit to 
#SBATCH --mem=250            # Memory pool for all cores 
#SBATCH -o enr.%A.%a.out       # File to which STDOUT will be written 
#SBATCH -e enr.%A.%a.err       # File to which STDERR will be written 

source new-modules.sh
module load bedtools2/2.25.0-fasrc01

python enrich_bed_sampler.py 5000 ${1##*/}_${2}_${3} $1 $2 $3 ${1##*/}_${2}_${3}_"${SLURM_ARRAY_TASK_ID}"
