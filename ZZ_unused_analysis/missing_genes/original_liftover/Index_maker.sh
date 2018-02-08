#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-1:00:00 #0-5:00:00 Runtime in days-hours:min:sec
#SBATCH --mem=100
#SBATCH -p serial_requeue #Partition to submit to
#SBATCH -o indexmaker%A.out
#SBTACH -e indexmaker%A.err

for FILE in $(ls *.fa)
do
	export FILE #this allows FILE to be passed to python with os.getenv
	python index_maker.py
done
