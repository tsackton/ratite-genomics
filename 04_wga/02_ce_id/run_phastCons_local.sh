#!/bin/bash

#define variables
VER=$1
EXPLEN=45
TARCOV=0.30

#set up
mkdir -p final_run_ver${VER}
cd final_run_ver$VER
cp ../final.ver${VER}.cons.mod cons.mod
cp ../final.ver${VER}.noncons.mod noncons.mod
mkdir -p logs
mkdir -p ELEMENTS SCORES
for INPUT in $(ls /n/regal/edwards_lab/ratites/wga/phast/final_mafs/split_all/*.ss)
do
	BASE=${INPUT##*/}
	SAMP=${BASE%%.ss}
	phastCons --expected-length=$EXPLEN --target-coverage=$TARCOV --most-conserved ELEMENTS/$SAMP.bed --score --msa-format SS $INPUT cons.mod,noncons.mod 1> ./SCORES/$SAMP.wig 2> ./logs/$SAMP.log &
done
cd ..
