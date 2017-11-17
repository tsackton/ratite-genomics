#!/bin/bash
#SBATCH -J pic_rem_dups
#SBATCH -n 1                     # Use 1 cores for the job
#SBATCH -t 0-04:00                 # Runtime in D-HH:MM
#SBATCH -p serial_requeue         # Partition to submit to
#SBATCH --mem=16200               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o pRd.%A.out  # File to which STDOUT will be written
#SBATCH -e pRd.%A.err  # File to which STDERR will be written

source new-modules.sh
module purge
module load picard/2.9.0-fasrc01

java -Xmx16g -jar ~/programs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT I=$1 O=no_dups_$1 M=remove_dup_metrics_$1
