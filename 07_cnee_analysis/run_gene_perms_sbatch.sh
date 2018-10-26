#!/usr/bin/bash

#SBATCH -t 480
#SBATCH --mem 8000
#SBATCH -p serial_requeue
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -e /dev/null
#SBATCH -o /dev/null
#SBATCH --array=1-100

source ~/default_env.sh
Rscript run_gene_perms.R galgal4 ${SLURM_ARRAY_TASK_ID} 8 100
Rscript run_gene_perms.R galgal5 ${SLURM_ARRAY_TASK_ID} 8 100