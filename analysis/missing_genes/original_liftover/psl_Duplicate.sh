#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-02:00:00 #0-0:05:00 Runtime in days-hours:min:sec
#SBATCH --mem 125
#SBATCH -p serial_requeue #Partition to submit to
#SBATCH -o psl_dup%A.out
#SBTACH -e psl_dup%A.err

for FILE in $(ls *.psl) 

do
        export FILE
        python psl_duplicate_parser.py
done
