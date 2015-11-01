#!/bin/bash

#these are the parameters that will vary with each iteration

export TARGETCOV=$1
export ESTLEN=$2

#get parts to use to estimate rho
cat part.* | chooseLines -k 100 - > random-sample.txt
split -a 3 -d -l 10 random-sample.txt sample. #make file parts
NUMFILES1=$(ls sample.* | wc -l)
NUMFILES=$(($NUMFILES1-1))

#submit estimate rho jobs and keep job id
rhoJobID=$(sbatch --array=0-$NUMFILES est_rho.sh | awk '{print $NF}')

#submit model averager job with dependency after rhoJobID
aveJobID=$(sbatch --dependency=afterany:$rhoJobID average_mods.sh | awk '{print $NF}')

#submit phastCons job
sbatch --dependency=afterok:$aveJobID run_phastCons.sh

