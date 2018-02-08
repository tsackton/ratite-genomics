#!/bin/bash
#SBATCH -J macs2
#SBATCH -n 1                     # Use 1 cores for the job
#SBATCH -t 0-12:00                 # Runtime in D-HH:MM
#SBATCH -p serial_requeue         # Partition to submit to
#SBATCH --mem=16000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o macs2.%A.out  # File to which STDOUT will be written
#SBATCH -e macs2.%A.err  # File to which STDERR will be written

source new-modules.sh
module purge
module load macs2/2.1.0.20140616-fasrc01

macs2 callpeak --name $1_macs2_pool_keepDup_q05 --treatment $1 $2 $3 --format BAMPE --q 0.05 --bdg --keep-dup all --gsize 1230260000
