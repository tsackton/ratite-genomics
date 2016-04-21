#!/bin/bash

#SBATCH -t 0-4:00
#SBATCH --mem 2000
#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=1-100

module load R/3.2.2-fasrc03
mkdir -p output/$1
Rscript --vanilla run_perms.R ${SLURM_ARRAY_TASK_ID} $1
