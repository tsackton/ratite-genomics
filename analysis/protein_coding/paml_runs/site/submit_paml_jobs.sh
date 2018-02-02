#!/bin/bash

module load paml/4.8-fasrc01

export HOGLIST=($(cat $1))
NUMHOG=${#HOGLIST[@]}
ZBNUMHOG=$(($NUMHOG - 1))
sbatch --array=0-$ZBNUMHOG ./paml_launcher_el.sh $1
