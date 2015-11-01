#!/bin/bash

#SBATCH -t 1-00:00
#SBATCH --mem 1500
#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=0-101
#SBATCH -J randAccel

RANDNUM=$1
BRST=$2
mkdir -p rand$RANDNUM
INDEX=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
for FILE in $(cat part.$INDEX)
do
	SAMP=${FILE%.ss}
	phyloP --method LRT --features LoweCNEEs.galGal4.bed.fixed --mode CONACC --branch ${BRST%?} $HOME/ratite_scratch/phast/neutMods/neut_ver1_final.named.mod $HOME/ratite_scratch/phast/phastCons/chunks/$SAMP.ss > ./rand$RANDNUM/$SAMP.features.out
done
