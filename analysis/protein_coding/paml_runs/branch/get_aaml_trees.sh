#!/bin/bash

HOG=$1
SUBDIRNUM=`expr ${HOG##HOG} % 100`
printf -v SUBDIR "%04d" $SUBDIRNUM
PAMLDIR=$SUBDIR/$HOG	
FINALDIR=$PAMLDIR/$HOG.codeml.aaml.ctl.out
OUTFILE=$FINALDIR/aaml.out
echo $OUTFILE
DATE=$(date +%F)
#tree
AATREE=$(grep -A4 "tree length" $OUTFILE | egrep "^[(]+[A-Za-z]")
IFS=$'\n'; arrDN=($AATREE); unset IFS
echo -e	"tree1\t$HOG\t${arrDN[0]}\n" >> aa_trees-$DATE.out
echo -e	"tree2\t$HOG\t${arrDN[1]}\n" >> aa_trees-$DATE.out
