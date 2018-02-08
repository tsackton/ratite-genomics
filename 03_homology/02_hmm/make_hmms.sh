#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 8000
#SBATCH --time 2-00:00
#SBATCH -e hmm_ogs.err
#SBATCH -o hmm_ogs.out
#SBATCH -J hmm_ogs

OUTPATH="hmms"
INPATH="alignments/forHMM"

# HMMER 3.1b2 (February 2015); http://hmmer.org/
# Copyright (C) 2015 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).

for FULLFILE in $(ls $INPATH/*.aligned)
do
	FILE=${FULLFILE##*/}
	ALN=${FILE%%.*}
	if [ ! -s "$OUTPATH/$ALN.hmm" ]
	then
		hmmbuild -n $ALN -o hmm_logs/$ALN.log $OUTPATH/$ALN.hmm $INPATH/$ALN.aligned 
	fi
done
