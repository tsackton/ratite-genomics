#!/bin/bash

#SBATCH -t 0-16:00
#SBATCH --mem 400
#SBATCH -p serial_requeue
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --array=0-111

INDEX=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
for SAMP in $(cat hmmpart.$INDEX)
do	
	if [ ! -e logs/$SAMP.done ]
	then
		#not done
		hmmsearch -E 1e-10 --cpu 8 -Z 19153963251 --tblout output/$SAMP.out --noali hmms/$SAMP.hmm seqs_for_hmm.fa > /dev/null	
		touch logs/$SAMP.done
	fi
done
