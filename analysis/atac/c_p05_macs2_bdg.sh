#!/bin/bash
#SBATCH -J mac2
#SBATCH -n 1                     # Use 1 cores for the job
#SBATCH -t 0-18:00                 # Runtime in D-HH:MM
#SBATCH -p serial_requeue         # Partition to submit to
#SBATCH --mem=20000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o macs2.%A.out  # File to which STDOUT will be written
#SBATCH -e macs2.%A.err  # File to which STDERR will be written

#$1 is name for file (and the bam file out of bowtie)

source new-modules.sh
module purge
module load macs2/2.1.0.20140616-fasrc01

macs2 callpeak --name $1_macs2_p05 --treatment $1 seq2_$1 --format BAMPE --p 0.05 --bdg --gsize 1230260000 --nomodel --extsize 200
