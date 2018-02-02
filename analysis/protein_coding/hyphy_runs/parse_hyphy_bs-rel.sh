#!/bin/bash

##HYPHY SCRIPT##
##all data in /n/holylfs/LABS/edwards_lab/tsackton/HyPhy/DATA

WORKDIR=$(pwd -P)
TESTCLASS=$1
DATE=$(date +%F)

for HOG in $(cat all_hogs)
do
    SUBDIRNUM=`expr ${HOG##HOG} % 100` 
    printf -v SUBDIR "%04d" $SUBDIRNUM
    RUNDIR=$SUBDIR/$HOG
    OUTPUT1=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/BS-REL/$RUNDIR/tree1/$HOG.aBSREL.OUT
    OUTPUT2=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/BS-REL/$RUNDIR/tree2/$HOG.aBSREL.OUT
    TREE1=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/BS-REL/$RUNDIR/tree1/$HOG.aBSREL.OUT.annotated.nwk
    TREE2=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/BS-REL/$RUNDIR/tree2/$HOG.aBSREL.OUT.annotated.nwk

    RESULTS1=$(python parse_bsrel_hog.py $HOG $TREE1 $OUTPUT1 $TESTCLASS)
    RESULTS2=$(python parse_bsrel_hog.py $HOG $TREE2 $OUTPUT2 $TESTCLASS)

    echo -e "$TESTCLASS\ttree1\t$RESULTS1" >> bsrel_res_parsed_${TESTCLASS}_${DATE}.txt
    echo -e "$TESTCLASS\ttree2\t$RESULTS2" >> bsrel_res_parsed_${TESTCLASS}_${DATE}.txt
done
