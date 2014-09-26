#!/bin/bash

#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 8000
#SBATCH -t 24:00:00
#SBATCH -J dedup
#SBATCH -o dedup_%j.out
#SBATCH -e dedup_%j.err


#sort frag and jump libraries, mark duplicates, and merge into a single library
#output _clean.bam is ready for SNP calling or other applications requiring deduplication

SPEC=$1

java -jar ~/sw/progs/picard-tools-1.121/SortSam.jar INPUT=${SPEC}_frag_mapped.bam OUTPUT=${SPEC}_sorted1.bam
java -jar ~/sw/progs/picard-tools-1.121/SortSam.jar INPUT=${SPEC}_jump_mapped.bam OUTPUT=${SPEC}_sorted2.bam
java -jar ~/sw/progs/picard-tools-1.121/MarkDuplicates.jar INPUT=${SPEC}_sorted1.bam OUTPUT=${SPEC}_dedup1.bam METRICS_FILE=${SPEC}_metrics1.txt REMOVE_DUPLICATES=true
java -jar ~/sw/progs/picard-tools-1.121/MarkDuplicates.jar INPUT=${SPEC}_sorted2.bam OUTPUT=${SPEC}_dedup2.bam METRICS_FILE=${SPEC}_metrics2.txt REMOVE_DUPLICATES=true
java -jar ~/sw/progs/picard-tools-1.121/MergeSamFiles.jar INPUT=${SPEC}_dedup1.bam INPUT=${SPEC}_dedup2.bam OUTPUT=${SPEC}_clean.bam

