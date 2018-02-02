#!/bin/bash
#SBATCH -J bdg2bw
#SBATCH -n 1                     # Use 1 cores for the job
#SBATCH -t 0-12:00                 # Runtime in D-HH:MM
#SBATCH -p serial_requeue         # Partition to submit to
#SBATCH --mem=70000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o bdg2bw.%A.out  # File to which STDOUT will be written
#SBATCH -e bdg2bw.%A.err  # File to which STDERR will be written

source new-modules.sh
module purge
module load bedtools2

for FILE in $(ls no*pool*bdg); do bedtools sort -i $FILE > sort_${FILE}; /n/home12/pgrayson/programs/bigWig/bedgraphToBigWig sort_${FILE} galGal.size ${FILE}.bw; done

