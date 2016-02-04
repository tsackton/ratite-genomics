#!/bin/bash

#define variables
VER=2
EXPLEN=45
TARCOV=0.30

#set up
mkdir -p final_run_ver${VER}_noratite
cd final_run_ver${VER}_noratite
#prune models for ratite-removed analysis
tree_doctor --prune rheAme,rhePen,strCam,droNov,casCas,aptHaa,aptOwe,aptRow ../final.ver${VER}.noncons.mod > noncons.mod
tree_doctor --prune rheAme,rhePen,strCam,droNov,casCas,aptHaa,aptOwe,aptRow ../final.ver${VER}.cons.mod > cons.mod
mkdir -p logs
mkdir -p ELEMENTS SCORES
for INPUT in $(ls /n/regal/edwards_lab/ratites/wga/phast/final_mafs/split_nr/*.ss)
do
	BASE=${INPUT##*/}
	SAMP=${BASE%%.ss}
	phastCons --expected-length=$EXPLEN --target-coverage=$TARCOV --most-conserved ELEMENTS/$SAMP.bed --score --msa-format SS $INPUT cons.mod,noncons.mod 1> ./SCORES/$SAMP.wig 2> ./logs/$SAMP.log &
done
