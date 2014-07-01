#!/bin/bash

#SBATCH -n 1
#SBATCH --mem 1000
#SBATCH -J extsra
#SBATCH -o extsra_%j.out
#SBATCH -e extsra_%j.err
#SBATCH -p serial_requeue
#SBATCH -t 9:00:00
#SBATCH --array 0-8

#extracts SRA files using a job array and a bash index
#move to the directory containing SRA files, and change the job array size as necessary

module load bio/sratoolkit.2.3.3-4

FILES=($(ls *.sra))
SRA=${FILES[$SLURM_ARRAY_TASK_ID]}
fastq-dump --split-files --gzip $SRA
rm $SRA
