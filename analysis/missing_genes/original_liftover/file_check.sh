#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-0:05:00 #0-5:00:00 Runtime in days-hours:min:sec
#SBATCH --mem=10
#SBATCH -p serial_requeue #Partition to submit to
#SBATCH -o Incomplete_Files/folder_%a.out

mkdir -p Incomplete_Files
for FILE in $(ls ../../Final/output/${SLURM_ARRAY_TASK_ID}/)
do
        if [ ! -f output_mfa/${SLURM_ARRAY_TASK_ID}/${FILE}.fasta ]
        then
                echo ${FILE}
        fi
done
