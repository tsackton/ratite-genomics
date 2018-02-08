#!/bin/bash

#SBATCH -t 0-16:00
#SBATCH --mem 2500
#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=0-101

VERSION=$1
mkdir -p logs
mkdir -p mods
INDEX=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
for FILE in $(cat part.$INDEX)
do
        SAMP=${FILE%.ss}
        STATUS=$(tail -n 1 ./logs/$SAMP.out)
        if [[ "$STATUS" != "Done." ]]
        then
		phastCons --expected-length=45 --target-coverage=0.3 --estimate-rho ./mods/$SAMP.rho --no-post-probs --msa-format SS --log ./logs/$SAMP.log ../chunks/$SAMP.ss ../../../neutMods/neut_ver${VERSION}_final.named.mod &> ./logs/$SAMP.out
	fi
done
