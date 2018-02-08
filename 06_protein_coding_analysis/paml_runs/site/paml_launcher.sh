#!/bin/bash

#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 3000
#SBATCH -t 2-00:00
#SBATCH -o ./logs/paml_%A_%a.out
#SBATCH -e ./logs/paml_%A_%a.err

HOGLIST=($(cat $1))
HOGRUN=${HOGLIST[$SLURM_ARRAY_TASK_ID]}
./run_paml.sh $HOGRUN
