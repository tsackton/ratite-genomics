#!/bin/bash

#SBATCH -n 16
#SBATCH --mem 1000
#SBATCH -J trimFastq
#SBATCH -o trim_%j.out
#SBATCH -e trim_%j.err
#SBATCH -p serial_requeue
#SBATCH -t 4:00:00

#trims Illumina adaptors from paired-end sequencing using Trimmomatic
#may need to adjust parameters for specific circumstances
#requires Trimmomatic 0.32 or greater
#requires all fastq files to be trimmed to use the same adaptor file and be in the same directory

#this is based on a local install of Trimmomatc; change $TRIMPATH as necessary

TRIMPATH=~/sw/progs/Trimmomatic-0.32/trimmomatic

#make adapter file; assumes that all fastq files to be trimmed use the same adaptor sequence and so doesn't 
if [ ! -s adapters.fa ]
then
	echo -e ">PrefixPE/1\nAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PrefixPE/2\nGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" > adapters.fa
fi

#specify the name of the first pair as the first argument on the command line and the base name for the output as the second argument

INFILE=$1
OUTFILE=$2

FILES=($(ls *.fastq.gz)) #all files to be trimmed need to be in one directory
FASTQ=${FILES[$SLURM_ARRAY_TASK_ID]}
java -jar $TRIMPATH PE -threads 16 -basein $INFILE -baseout $OUTFILE ILLUMINACLIP:adapters.fa:2:30:10
