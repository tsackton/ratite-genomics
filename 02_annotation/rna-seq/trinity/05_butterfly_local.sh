#!/bin/bash
#SBATCH -p general # Partition to submit to
#SBATCH –n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH –t 1-0:00 # Runtime in days-hours:minutes
#SBATCH --mem 4000 # Memory in MB
#SBATCH –J trinity_b # job name
#SBATCH -o A1_trinity_b.out # File to which standard out will be written
#SBATCH –e A1_trinity_b.err # File to which standard err will be written

#get species to process
SP=$1
LAUNCHPATH=$(pwd -P)

# rename and compress again
LEFT=$(ls -m $LAUNCHPATH/$SP/trimmed/*.R1.*.fastq.gz | perl -p -i -e 's/\s+//')
RIGHT=$(ls -m $LAUNCHPATH/$SP/trimmed/*.R2.*.fastq.gz | perl -p -i -e 's/\s+//')

Trinity --seqType fq \
 --left $LEFT \
 --right $RIGHT \
 --output $SP/trinity \
 --max_memory 10G
