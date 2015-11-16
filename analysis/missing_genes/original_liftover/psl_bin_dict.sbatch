#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1-23:00:00 #0-5:00:00 Runtime in days-hours:min:sec
#SBATCH --mem=1500
#SBATCH -p serial_requeue #Partition to submit to
#SBATCH -o dic/%A.out
#SBATCH -e dic/%A.err

python psl_bin_dic_parser.py $1 $2
