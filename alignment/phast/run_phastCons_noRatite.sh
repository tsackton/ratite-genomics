#!/bin/bash

#define variables
VER=2
EXPLEN=45
TARCOV=0.30

#set up
mkdir -p logs
mkdir -p ELEMENTS SCORES
for INPUT in $(ls /n/regal/edwards_lab/ratites/phast/final_nor_mafs/*.ss)
do
	BASE=${INPUT##*/}
	SAMP=${BASE%%.ss}
	phastCons --expected-length=$EXPLEN --target-coverage=$TARCOV --most-conserved ELEMENTS/$SAMP.bed --score --msa-format SS $INPUT cons.mod,noncons.mod 1> ./SCORES/$SAMP.wig 2> ./logs/$SAMP.log &
done
