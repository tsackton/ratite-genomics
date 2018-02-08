#!/bin/bash
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-2:30 # Runtime in days-hours:minutes
#SBATCH --mem 500 # Memory in MB
#SBATCH -J A1_fastqc2 #job name
#SBATCH -o A1_fastqc2.out # File to which standard out will be written
#SBATCH -e A1_fastqc2.err # File to which standard err will be written

#get work dir
WORKDIR=${1:?} #directory to process, exit if not set

#get files
for FQ in $(ls $WORKDIR/trimmed/*.fastq.gz);
do
	fastqc -o $WORKDIR/fastqc_reports $FQ
done

