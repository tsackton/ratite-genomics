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
LIB=$2
RUN=$3

READ1="../data/trimmmed/${GENOME}_${RUN}_${LIB}_1P.fastq.gz"
READ2="../data/trimmmed/${GENOME}_${RUN}_${LIB}_2P.fastq.gz"
OUTPUT="$GENOME.$LIB.$RUN"

bwa mem -t 16 -M -R '@RG\tID:'"$GENOME"'\tSM:'"$GENOME"'\tPL:Illumina_'"$RUN"'\tLB:'"$LIB"'\tPU:'"$LIB.$RUN" "${GENOME}1" $READ1 $READ2 | samtools -bT ${GENOME}1.fa.gz - > ${OUTPUT}.bam
scontrol update JobId=$SLURM_JOBID NumNodes=1
java -Xmx2g -jar ~/sw/progs/picard-tools-1.121/SortSam.jar TMP_DIR=/scratch INPUT=${OUTPUT}.bam OUTPUT=${OUTPUT}.sorted.bam SORT_ORDER=coordinate
if [ $? -eq 0 ]
then
	rm ${OUTPUT}.bam
fi
java -Xmx2g -jar ~/sw/progs/picard-tools-1.121/MarkDuplicates.jar TMP_DIR=/scratch INPUT=${OUTPUT}.sorted.bam OUTPUT=${OUTPUT}.dedup.bam METRICS_FILE=${OUTPUT}.metrics.txt REMOVE_DUPLICATES=true
