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
READ1=$2
READ2=$3

bwa mem -t 16 -M $GENOME $READ1 $READ2 > ${GENOME}_frag_mapped.sam