#!/bin/bash

BASE=/n/regal/edwards_lab/ratites/maker2
RNASEQ=$BASE/evidence/rnaseq
PROT=$BASE/evidence/prot
SNAPDIR=$BASE/evidence/snap
GENOMEDIR=$BASE/inputs
RUNDIR=$BASE/annotation

mkdir -p $RNASEQ
mkdir -p $PROT
mkdir -p $SNAPDIR
mkdir -p $GENOMEDIR

cp $HOME/ratite_store/trinity/final/*final.fasta $RNASEQ
mv $BASE/rnaseq/*.gff $RNASEQ
cp $HOME/ratite_store/ratite-genomics/annotation/maker/run1/get_proteins.sh $PROT
cp $HOME/ratite_store/genomes/ASSEMBLY_FREEZE/*.fa.gz $GENOMEDIR
cp $HOME/ratite_scratch/training/*/*.hmm $SNAPDIR

#generate empty control files
maker -CTL

#bops and exe just use defaults so can be copied to species directories immediately

for SPECIES in aptHaa aptOwe aptRow droNov casCas notPer cryCin eudEle rheAme rhePen;
do
	mkdir -p $RUNDIR/$SPECIES
	cp maker_bopts.ctl $RUNDIR/$SPECIES
	cp maker_exe.ctl $RUNDIR/$SPECIES
	cp maker_opts.ctl maker_opts_$SPECIES.ctl
done

rm maker_bopts.ctl
rm maker_exe.ctl
rm maker_opts.ctl 

#the opts files for each species need to be manually edited unfortunately, although changes will be verified with diff






