#!/bin/bash

#SBATCH -p serial_requeue
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH -J tophat
#SBATCH -e tophat_%j.err
#SBATCH -o tophat_%j.out
#SBATCH --mem 16000
#SBATCH -A edwards_lab

source new-modules.sh
module load bowtie2/2.2.2-fasrc01
module load tophat/2.0.13-fasrc01

TARGET=$1
LEFT=$2
RIGHT=$3

##setup options, as they are different for kiwi and notper/dronov

TMP=${LEFT##reads/}
SP=${TMP%%/*}

OUTDIR="$TARGET.$SP"

#echo info to standard out
echo "Running tophat with reads from $SP aligning to $TARGET in output directory $OUTDIR"

if [ -d $OUTDIR ]
then
	tophat -R $OUTDIR
else
	if [ $SP == "kiwi" ]
	then
		tophat -N 8 --read-gap-length 3 --read-edit-dist 9 -i 20 --b2-very-sensitive --library-type fr-unstranded --no-coverage-search -p 8 -o $OUTDIR $TARGET $LEFT $RIGHT
	else
		tophat -N 8 --read-gap-length 3 --read-edit-dist 9 -i 20 --b2-very-sensitive --library-type fr-secondstrand --no-coverage-search -p 8 -o $OUTDIR $TARGET $LEFT $RIGHT
	fi
fi
