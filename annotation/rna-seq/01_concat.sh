#!/bin/bash
#SBATCH -p general # Partition to submit to
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-4:00 # Runtime in days-hours:minutes
#SBATCH --mem 1000 # Memory in MB
#SBATCH -J A01_fastqc #job name
#SBATCH -o logs/A01_fastqc_%j.out # File to which standard out will be written
#SBATCH -e logs/A01_fastqc_%j.err # File to which standard err will be written

#concatenate all left and right reads in a directory then QC
#first command line argument is the directory containing the files to process
#attempts to be smart about figuring out which is left and which is right read
#but you need to make sure that your files are named in a way the program can figure out

#command line variables
SEQDIR=${1:?} #directory to process, exit if not set
WORKDIR=${2:?} #working directory
SPNAME=${3:$WORKDIR} #species name, required for output; also assumed to be working dir

#set up dirs
mkdir -p $WORKDIR/raw
mkdir -p $WORKDIR/fastqc_reports

# concat first
cat $SEQDIR/*1.fastq.gz > $WORKDIR/raw/$SPNAME.R1.fastq.gz
cat $SEQDIR/*2.fastq.gz > $WORKDIR/raw/$SPNAME.R2.fastq.gz

# and now QC
fastqc -o $WORKDIR/fastqc_reports $WORKDIR/raw/$SPNAME.R1.fastq.gz
fastqc -o $WORKDIR/fastqc_reports $WORKDIR/raw/$SPNAME.R2.fastq.gz
