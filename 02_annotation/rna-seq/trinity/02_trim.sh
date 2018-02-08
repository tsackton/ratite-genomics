#!/bin/bash
#SBATCH -p general # Partition to submit to
#SBATCH -n 8 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-12:00 # Runtime in days-hours:minutes
#SBATCH --mem 12000 # Memory in MB
#SBATCH -J A02_trim # job name
#SBATCH -o logs/A02_trim_%j.out # File to which standard out will be written
#SBATCH -e logs/A02_trim_%j.err # File to which standard err will be written

module load legacy
module load centos6/Trimmomatic-0.32 

#command line variables
SEQDIR=${1:?} #directory to process, exit if not set; assumes that raw sequence is in raw subdir
SPNAME=${2:?} #species name, required for output; exit if not set

#set up dirs
mkdir -p $SEQDIR/trimmed

cat $TRIMMOMATIC/adapters/TruSeq2-PE.fa $TRIMMOMATIC/adapters/TruSeq3-PE-2.fa > adapters.fa
java -jar $TRIMMOMATIC/trimmomatic-0.32.jar PE \
 -threads 8 \
 $SEQDIR/raw/$SPNAME.R1.fastq.gz $SEQDIR/raw/$SPNAME.R2.fastq.gz \
 $SEQDIR/trimmed/$SPNAME.R1.pair.fastq.gz $SEQDIR/trimmed/$SPNAME.R1.single.fastq.gz \
 $SEQDIR/trimmed/$SPNAME.R2.pair.fastq.gz $SEQDIR/trimmed/$SPNAME.R2.single.fastq.gz \
 ILLUMINACLIP:adapters.fa:2:30:10:1:true \
 LEADING:3 TRAILING:3 \
 SLIDINGWINDOW:4:10 MINLEN:25
