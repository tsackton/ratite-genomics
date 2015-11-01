#!/bin/bash

#SBATCH -t 0-16:00
#SBATCH --mem 2500
#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1

mkdir -p logs
mkdir -p mods
INDEX=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
for FILE in $(cat sample.$INDEX)
do
        SAMP=${FILE%.ss}
        STATUS=$(tail -n 1 ./logs/$SAMP.out)
        if [[ "$STATUS" != "Done." ]]
        then
		phastCons --expected-length=$ESTLEN --target-coverage=$TARGETCOV --estimate-rho ./mods/$SAMP.rho --no-post-probs --msa-format SS --log ./logs/$SAMP.log ../chunks/$SAMP.ss ../../neutMods/neut_ver1_final.named.mod &> ./logs/$SAMP.out
	fi
done
