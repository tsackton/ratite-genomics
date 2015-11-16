#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-00:30:00 #0-0:05:00 Runtime in days-hours:min:sec
#SBATCH --mem 20
#SBATCH -p serial_requeue #Partition to submit to
#SBATCH -o psl_checker%A.out
#SBTACH -e psl_checker%A.err

for FILE in $(ls *.psl) 

do
	export FILE
	python psl_tester.py
done
