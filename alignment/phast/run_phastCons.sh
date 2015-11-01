#!/bin/bash

#SBATCH -t 1-00:00
#SBATCH --mem 4000
#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=0-101

mkdir -p logs
mkdir -p ELEMENTS SCORES
INDEX=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
for FILE in $(cat part.$INDEX)
do
        SAMP=${FILE%.ss}
	STATUS=$(tail -n 1 logs/$SAMP.log)
        if [[ "$STATUS" != "Done." ]]
        then
		phastCons --expected-length=$ESTLEN --target-coverage=$TARGETCOV --most-conserved ELEMENTS/$SAMP.bed --score --msa-format SS ../chunks/$SAMP.ss ave.cons.mod,ave.noncons.mod 1> ./SCORES/$SAMP.wig 2> ./logs/$SAMP.log
	fi
done
