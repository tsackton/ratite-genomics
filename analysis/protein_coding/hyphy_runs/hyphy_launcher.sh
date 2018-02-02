#!/bin/bash

#SBATCH -p general
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --constraint="holyib"
#SBATCH -A informatics
#SBATCH --mem 8000
#SBATCH -t 7-0:00
#SBATCH -o ./logs/paml_%A_%a.out
#SBATCH -e ./logs/paml_%A_%a.err

HOGLIST=($(cat $1))
HOGRUN=${HOGLIST[$SLURM_ARRAY_TASK_ID]}
./run_hyphy_bs-rel.sh $HOGRUN
