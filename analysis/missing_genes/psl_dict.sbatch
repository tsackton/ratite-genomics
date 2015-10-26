#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-12:00:00 #0-5:00:00 Runtime in days-hours:min:sec
#SBATCH --mem=1000
#SBATCH -p serial_requeue #Partition to submit to
#SBATCH -o dic/%a.out
#SBATCH -e dic/%a.err

python psl_dictionary_parser.py $1
