#!/bin/bash

#SBATCH -p tsackton,serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 1000
#SBATCH -t 0-18:00
#SBATCH -o ./logs/paml_%A_%a.out
#SBATCH -e ./logs/paml_%A_%a.err

HOGLIST=($(cat rerun.round2))
HOGRUN=${HOGLIST[$SLURM_ARRAY_TASK_ID]}
./run_paml.sh $HOGRUN
