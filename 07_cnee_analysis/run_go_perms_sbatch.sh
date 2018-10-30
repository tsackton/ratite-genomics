#!/usr/bin/bash

#SBATCH -t 800
#SBATCH --mem 40000
#SBATCH -p serial_requeue
#SBATCH -n 5
#SBATCH -N 1
#SBATCH --array=1-100

source ~/default_env.sh
Rscript run_enrichment_perms.R galgal4 ${SLURM_ARRAY_TASK_ID} 5 50
Rscript run_enrichment_perms.R galgal5 ${SLURM_ARRAY_TASK_ID} 5 50
