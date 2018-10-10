#!/usr/bin/bash

#SBATCH -t 120
#SBATCH --mem 4000
#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -e /dev/null
#SBATCH -o /dev/null

source ~/default_env.sh
Rscript --vanilla $1 $2 $3