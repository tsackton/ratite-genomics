#!/usr/bin/bash

#SBATCH -t 480
#SBATCH --mem 40000
#SBATCH -p serial_requeue
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --array=1-100

source ~/default_env.sh
Rscript run_convergence_perms.R ${SLURM_ARRAY_TASK_ID} 4 50
