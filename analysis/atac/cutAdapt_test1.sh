#!/bin/bash
#SBATCH -J cutadapt
#SBATCH -n 1                     # Use 1 cores for the job
#SBATCH -t 0-9:00                 # Runtime in D-HH:MM
#SBATCH -p serial_requeue         # Partition to submit to
#SBATCH --mem=100               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o cutad.%A.out  # File to which STDOUT will be written
#SBATCH -e cutad.%A.err  # File to which STDERR will be written

source new-modules.sh
module purge
module load cutadapt/1.8.1-fasrc01
module load python/2.7.6-fasrc01

# $1 = R1 reads
# $2 = R2 reads

cutadapt -b CTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAGGGGG -b CTGTCT -B CTGTCTCTTATACACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAGGG -o trim_test1_$1 -p trim_test1_$2  $1 $2
