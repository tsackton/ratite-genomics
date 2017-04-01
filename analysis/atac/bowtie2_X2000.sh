#!/bin/bash

# Call this script with 2 fastq file2 for reads 1 and 2.
#
# e.g. sbatch [thisscript.sh](http://thisscript.sh) myfile.R1.fastq  myfile.R2.fastq
#
# myfile.R1.fastq is referenced by the variable $1
# myfile.R2.fastq is referenced by the variable $2

#SBATCH -J ATAC_Bowtie2 
#SBATCH -N 1                      # Ensure that all cores are on one machine
#SBATCH -n 32                # Use n cores for one job 
#SBATCH -t 0-12:00                # Runtime in D-HH:MM 
#SBATCH -p serial_requeue                # Partition to submit to 
#SBATCH --mem=32000            # Memory pool for all cores 
#SBATCH -o bt2ATAC.%A.out       # File to which STDOUT will be written 
#SBATCH -e bt2ATAC.%A.err       # File to which STDERR will be written 

source new-modules.sh
module purge
module load bowtie2/2.2.4-fasrc01
module load samtools/0.1.19-fasrc01

bowtie2 -x $3 -1 $1 -2 $2 -X 2000 -p 32 | samtools view -b -S - |samtools sort - $1

samtools index $1.bam
