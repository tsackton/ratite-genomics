#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem 8000
#SBATCH --time 2-00:00
#SBATCH -e aln_ogs.err
#SBATCH -o aln_ogs.out
#SBATCH -J aln_ogs

INPATH="OutputFolder/HOGFasta"
OUTPATH="alignments/forHMM"

#
#------------------------------------------------------------------------------
#  MAFFT v7.221 (2014/04/16)
#  http://mafft.cbrc.jp/alignment/software/
#  MBE 30:772-780 (2013), NAR 30:3059-3066 (2002)
#------------------------------------------------------------------------------
#

for FULLFILE in $(ls $INPATH/*.fa)
do
	FILE=${FULLFILE##*/}
	ALN=${FILE%%.*}
	if [ ! -s "$OUTPATH/$ALN.aligned" ]
	then
		mafft --globalpair --maxiterate 1000 --thread 8 $FULLFILE > $OUTPATH/$ALN.aligned
	fi
done
