#!/bin/bash

#SBATCH -p general
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mem 2500
#SBATCH -J removePhiX
#SBATCH -o removePhix_%j.out
#SBATCH -e removePhix_%j.err
#SBATCH -t 4:00:00
 
$HOME/sw/progs/bbmap/bbduk.sh -Xmx2g in=GSK_R1_220demultiplex2.fastq in2=GSK_R2_220demultiplex2.fastq out=GSK_R1_220demultiplex2_nophix.fastq.gz out2=GSK_R2_220demultiplex2_nophix.fastq.gz \
   stats=gsk_stats.txt showspeed=f qin=33 qout=33 threads=8 k=31 ref=phix.fa edist=1

