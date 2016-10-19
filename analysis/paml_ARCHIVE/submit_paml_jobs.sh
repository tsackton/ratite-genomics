#!/bin/bash

export HOGLIST=($(cat $1))
NUMHOG=${#HOGLIST[@]}
ZBNUMHOG=$(($NUMHOG - 1))
sbatch --array=0-$ZBNUMHOG ./paml_launcher.sh $1
