#!/bin/bash

#SBATCH -t 0-08:00
#SBATCH --mem 1500
#SBATCH -p serial_requeue
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=0-101

mkdir -p ratite
mkdir -p tinamou
INDEX=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
for FILE in $(cat part.$INDEX)
do
	SAMP=${FILE%.ss}
	phyloP --method LRT --features LoweCNEEs.galGal4.bed.fixed --mode CONACC --branch rheAme,rhePen,rheAme-rhePen,strCam,casCas,droNov,casCas-droNov,aptHaa-casCas,aptHaa-aptOwe,aptRow,aptHaa,aptOwe,aptHaa-aptRow $HOME/ratite_scratch/phast/neutMods/neut_ver1_final.named.mod $HOME/ratite_scratch/phast/phastCons/galGal/chunks/$SAMP.ss > ./ratite/$SAMP.features.ratite.out
	phyloP --method LRT --features LoweCNEEs.galGal4.bed.fixed --mode CONACC --branch cryCin,tinGut,cryCin-tinGut,eudEle,notPer,eudEle-notPer,cryCin-eudEle $HOME/ratite_scratch/phast/neutMods/neut_ver1_final.named.mod $HOME/ratite_scratch/phast/phastCons/galGal/chunks/$SAMP.ss > ./tinamou/$SAMP.features.tinamou.out
done
