#!/bin/bash

#SBATCH -n 16
#SBATCH --mem 16000
#SBATCH -J trimFastq
#SBATCH -o trim_%j.out
#SBATCH -e trim_%j.err
#SBATCH -p general
#SBATCH -t 4:00:00

#trims Illumina adaptors from paired-end sequencing using Trimmomatic
#may need to adjust parameters for specific circumstances
#requires Trimmomatic 0.32 or greater
#requires all fastq files to be trimmed to use the same adaptor file and be in the same directory

#this is based on a local install of Trimmomatc; change $TRIMPATH as necessary

TRIMPATH=~/sw/progs/Trimmomatic-0.32/trimmomatic-0.32.jar

#make adapter file; assumes that all fastq files to be trimmed use the same adaptor sequence
#note that this current file merges nextera and other adapters which should be okay but may not be ideal for final processing
#nextera sequences are derived from files distributed with Trimmomatic but look correct
if [ ! -s adapters.fa ]
then
	echo -e ">PrefixPE/1\nAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PrefixPE/2\nGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" > adapters.fa
	echo -e ">PrefixNX/1\nAGATGTGTATAAGAGACAG\n>PrefixNX/2\nAGATGTGTATAAGAGACAG" >> adapters.fa
	echo -e ">Trans1\nTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG\n>Trans1_rc\nCTGTCTCTTATACACATCTGACGCTGCCGACGA\n>Trans2\nGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG\n>Trans2_rc\nCTGTCTCTTATACACATCTCCGAGCCCACGAGAC\n" >> adapters.fa
fi

#specify the name of the first pair as the first argument on the command line and the base name for the output as the second argument

INFILE=$1
OUTFILE=$2

FILES=($(ls *.fastq.gz)) #all files to be trimmed need to be in one directory
FASTQ=${FILES[$SLURM_ARRAY_TASK_ID]}
java -jar $TRIMPATH PE -threads 16 -basein $INFILE -baseout $OUTFILE ILLUMINACLIP:adapters.fa:2:30:10:1:true