#!/bin/bash
#SBATCH -p general # Partition to submit to
#SBATCH –n 16 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH –t 3-0:00 # Runtime in days-hours:minutes
#SBATCH --mem 155000 # Memory in MB
#SBATCH –J norm_ic # job name
#SBATCH -o norm_ic_%j.out # File to which standard out will be written
#SBATCH –e norm_ic_%j.err # File to which standard err will be written

#this script normalizes the properly paired reads only, does not deal with single-end reads
#source new-modules.sh
#module load trinityrnaseq

#get species to process
SP=$1
LAUNCHPATH=$(pwd -P)

# rename and compress again
LEFT=$(ls -m $LAUNCHPATH/$SP/trimmed/*.R1.*.fastq.gz | perl -p -i -e 's/\s+//')
RIGHT=$(ls -m $LAUNCHPATH/$SP/trimmed/*.R2.*.fastq.gz | perl -p -i -e 's/\s+//')

# normalize and run Inchworm + Chryslalis
Trinity \
 --seqType fq \
 --max_memory 500G \
 --left $LEFT \
 --right $RIGHT \
 --normalize_reads \
 --output $SP/trinity \
 --CPU 8 \
 --verbose \
 --no_distributed_trinity_exec \
 --min_kmer_cov 1 \
 --group_pairs_distance 800
