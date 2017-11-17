#!/bin/bash
#SBATCH -J mito_rem
#SBATCH -n 1                     # Use 1 cores for the job
#SBATCH -t 0-04:00                 # Runtime in D-HH:MM
#SBATCH -p serial_requeue         # Partition to submit to
#SBATCH --mem=250               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o mitoR.%A.out  # File to which STDOUT will be written
#SBATCH -e mitoR.%A.err  # File to which STDERR will be written

source new-modules.sh
module purge
module load samtools
module load ATAC-seq

samtools view -h $1 | removeChrom - - NC_001323.1 | samtools view -b - > noMito_$1
