#!/bin/bash

#SBATCH -t 0-24:00
#SBATCH --mem 4000
#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=100-102

module load R
Rscript --vanilla run_perms.R ${SLURM_ARRAY_TASK_ID}
