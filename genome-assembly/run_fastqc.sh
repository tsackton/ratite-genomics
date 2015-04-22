#!/bin/bash

#SBATCH -p serial_requeue
#SBATCH --mem 1000
#SBATCH -t 600
#SBATCH -n 1
#SBATCH --array=0-21
#SBATCH -J fastqc
#SBATCH -o fastqc_%j.out
#SBATCH -e fastqc_%j.err

module load centos6/fastqc-0.10.1

#need to adjust array to appropriate number
#assumes all files have been processed by trim script and are in one directory

FILES=($(ls *P.fastq.gz)) #gets the paired reads after trimming - creates list of files
FQ=${FILES[$SLURM_ARRAY_TASK_ID]} #selects file to run in each array job

fastqc --noextract -t 1 -k 5 -q $FQ

