#!/usr/bin/bash

#SBATCH -t 480
#SBATCH --mem 40000
#SBATCH -p serial_requeue,bos-info
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --array=1-100

source ~/default_env.sh
DATAPATH=/n/holylfs/LABS/edwards_lab/tsackton/RATITE_PAPER_DATA_FREEZE/DRYAD/07_cnees/processed
Rscript run_gene_perms.R galgal4 ${SLURM_ARRAY_TASK_ID} 8 100 $DATAPATH
Rscript run_gene_perms.R galgal5 ${SLURM_ARRAY_TASK_ID} 8 100 $DATAPATH
