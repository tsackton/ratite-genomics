#!/bin/bash

#SBATCH -p general
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --mem 8000
#SBATCH -t 36:00:00
#SBATCH -J bwa_mem
#SBATCH -o bwa_mem_%j.out
#SBATCH -e bwa_mem_%j.err

GENOME=$1
OUTPUT=$2
READ1=$3
READ2=$4

bwa mem -t 16 -M $GENOME $READ1 $READ2 | samtools view -bT $GENOME.fa.gz -F 0x4 -f 0x2 - > ${OUTPUT}.bam
