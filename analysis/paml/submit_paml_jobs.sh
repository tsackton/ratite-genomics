#!/bin/bash

export HOGLIST=($(cat rerun.round2))
NUMHOG=${#HOGLIST[@]}
ZBNUMHOG=$(($NUMHOG - 1))
sbatch --array=0-$ZBNUMHOG ./paml_launcher.sh
