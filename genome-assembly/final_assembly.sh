#!/bin/bash

#Final assembly script

#Define species to run
SPEC=$1

#Make sure we are in the right directory
cd /scratch/tsackton/allpaths_runs

#Prep paths
export PATH=$PATH:$HOME/sw/progs/allpathslg-50191
export LD_LIBRARY_PATH=/n/sw/graphviz-2.22.2/lib:/n/sw/fasrcsw/apps/Core/gcc/4.8.2-fasrc01/lib64:/n/sw/fasrcsw/apps/Core/gcc/4.8.2-fasrc01/lib64:/n/sw/fasrcsw/apps/Core/gcc/4.8.2-fasrc01/lib:/n/sw/centos6/mpc-1.0.1/lib:/n/sw/centos6/mpfr-3.1.2/lib:/n/sw/centos6/gmp-5.1.1/lib:/n/sw/fasrcsw/apps/Core/gmp/5.1.3-fasrc01/lib64:/lsf/7.0/linux2.6-glibc2.3-x86_64/lib

#Make directory structure
mkdir -p $SPEC/RAWSEQ/20141223

#Prepare inputs
if [[ -s $SPEC_groups.csv && -s $SPEC_libs.csv ]] 
then
	PrepareAllPathsInputs.pl DATA_DIR=/scratch/tsackton/allpaths_runs/$SPEC/RAWSEQ PLOIDY=2 IN_GROUPS_CSV=$SPEC_groups.csv IN_LIBS_CSV=$SPEC_libs.csv 1> $SPEC_prep.log 2> $SPEC_prep.err
else
	echo "Groups or Libs CSV is missing!"
	exit 1
fi

#Run ALLPATHS
RunAllPathsLG PRE=/scratch/tsackton/allpaths_runs REFERENCE_NAME=$SPEC DATA_SUBDIR=RAWSEQ RUN=20141223 HAPLOIDIFY=TRUE THREADS=30 OVERWRITE=TRUE 1> $SPEC_20141223.log

#Clean up
mkdir /scratch/tsackton/assembles/$SPEC
mv $SPEC_* /scratch/tsackton/assembles/$SPEC
mv run_$SPEC_haplo.sh /scratch/tsackton/assembles/$SPEC
mv $SPEC/make_log/ /scratch/tsackton/assembles/$SPEC
mv $SPEC/RAWSEQ/20141223/*.kspec /scratch/tsackton/assembles/$SPEC
mv $SPEC/RAWSEQ/20141207/ASSEMBLIES/test/final.{assembly,contigs}.{fasta,fastg,efasta} /scratch/tsackton/assembles/$SPEC
mv $SPEC/RAWSEQ/20141207/ASSEMBLIES/test/*.report /scratch/tsackton/assembles/$SPEC
rm -r $SPEC
