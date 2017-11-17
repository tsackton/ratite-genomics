#!/bin/bash
#SBATCH -J NGmerge
#SBATCH -n 4                     # Use 1 cores for the job
#SBATCH -t 0-04:00                 # Runtime in D-HH:MM
#SBATCH -p serial_requeue         # Partition to submit to
#SBATCH --mem=100               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o NGmerge.%A.out  # File to which STDOUT will be written
#SBATCH -e NGmerge.%A.err  # File to which STDERR will be written

source new-modules.sh
module purge
module load ATAC-seq

NGmerge -a -n 4 -e 20 -v -1 $1.R1.fastq.gz -2 $1.R2.fastq.gz -o 20_$1

