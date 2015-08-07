#!/bin/bash

#SBATCH -n 1
#SBATCH --mem 1000
#SBATCH -J extsra
#SBATCH -o logs/A00_extsra_%j.out
#SBATCH -e logs/A00_extsra_%j.err
#SBATCH -p general
#SBATCH -t 10:00:00
#SBATCH --array 0-30 #large array size; array jobs larger than number of files will exit

#extracts SRA files using a job array and a bash index
#first command line argument is the directory containing the SRA files to process
#second command line argument is 1 to remove SRA files after they are extracted and 0 to keep them

module load legacy
module load bio/sratoolkit.2.3.3-4

#command line variables
SRADIR=${1:?} #directory to process, exit if not set
REMOVE=${2:-0} #0 to keep SRA files, 1 to remove, default is keep if this parameter is not set

#move to SRA directory
cd $SRADIR

#create bash array with all SRA files in direcory
FILES=($(ls *.sra))

#get size of array
NUMSRA=${#FILES[@]}

#check to make sure that only jobs with valid files run
if [[ $SLURM_ARRAY_TASK_ID -ge $NUMSRA ]]
then
	echo "Array out of bounds"
	exit
fi

#actual processing of sra files
SRA=${FILES[$SLURM_ARRAY_TASK_ID]}
fastq-dump --split-files --gzip $SRA
if [[$REMOVE -gt 0 ]]
then
	rm $SRA
fi

#done
