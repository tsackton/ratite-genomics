#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-0:30:00 #0-5:00:00 Runtime in days-hours:min:sec
#SBATCH --mem=2000
#SBATCH -p serial_requeue #Partition to submit to
#SBATCH -o mfr/%A.out
#SBATCH -e mfr/%A.err

python stop_number_3frame.py $1
