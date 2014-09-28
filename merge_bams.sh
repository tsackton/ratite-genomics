#!/bin/bash

#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 8000
#SBATCH -t 24:00:00
#SBATCH -J merge
#SBATCH -o mergebams_%j.out
#SBATCH -e mergebams_%j.err

#merge output for each species

SPEC=$1
RAW=$(ls $SPEC*.sorted.bam)
CLEAN=$(ls $SPEC*.dedup.bam)

samtools merge $SPEC.raw.bam $RAW 
if [ $? -eq 0 ]
then
	rm $SPEC*.sorted.bam
fi
samtools merge -u - $CLEAN | samtools view -F 0x4 -f 0x2 -b > $SPEC.tmp.bam
java -Xmx2g -jar ~/sw/progs/picard-tools-1.121/MarkDuplicates.jar TMP_DIR=/scratch INPUT=$SPEC.tmp.bam OUTPUT=${SPEC}.clean.bam METRICS_FILE=${SPEC}.clean.metrics.txt REMOVE_DUPLICATES=true
if [ $? -eq 0 ]
then
	rm ${SPEC}.tmp.bam
	rm $SPEC*.dedup.bam
fi