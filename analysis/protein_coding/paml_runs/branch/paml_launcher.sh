#!/bin/bash

#SBATCH -p shared,general
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 3000
#SBATCH -t 1-00:00
#SBATCH -o ./logs/aaml_%A_%a.out
#SBATCH -e ./logs/aaml_%A_%a.err

HOGLIST=($(cat $1))
HOGRUN=${HOGLIST[$SLURM_ARRAY_TASK_ID]}
./run_paml_aa.sh $HOGRUN
