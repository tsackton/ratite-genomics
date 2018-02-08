#!/bin/bash

##HYPHY SCRIPT##
##all data in /n/holylfs/LABS/edwards_lab/tsackton/HyPhy/DATA

WORKDIR=$(pwd -P)
HOG=$1
SUBDIRNUM=`expr ${HOG##HOG} % 100`
printf -v SUBDIR "%04d" $SUBDIRNUM
RUNDIR=$SUBDIR/$HOG
mkdir -p $RUNDIR
cd $RUNDIR

#set up run
mkdir -p tree1
mkdir -p tree2
ALIGN=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/DATA/$RUNDIR/$HOG.phy
TREE1=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/DATA/$RUNDIR/$HOG.tree1.nwk
TREE2=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/DATA/$RUNDIR/$HOG.tree2.nwk
#cleanup trees and alignment
nw_topology $TREE1 | perl -p -e 's/[.-]/_/g' > tree1/$HOG.tree1.nwk
nw_topology $TREE2 | perl -p -e 's/[.-]/_/g' > tree2/$HOG.tree2.nwk
cat $ALIGN | perl -p -e 's/[.-]/_/g' > $HOG.phy
#make paths for input
ALIGNIN=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/BS-REL/$RUNDIR/$HOG.phy
TR1IN=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/BS-REL/$RUNDIR/tree1/$HOG.tree1.nwk
TR2IN=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/BS-REL/$RUNDIR/tree2/$HOG.tree2.nwk
OUTPUT1=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/BS-REL/$RUNDIR/tree1/$HOG.aBSREL.OUT
OUTPUT2=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/BS-REL/$RUNDIR/tree2/$HOG.aBSREL.OUT
echo -e "1\n1\n2\n$ALIGNIN\n$TR1IN\n2\nd\n$OUTPUT1\n" > tree1/input
echo -e "1\n1\n2\n$ALIGNIN\n$TR2IN\n2\nd\n$OUTPUT2\n" > tree2/input

echo "$(pwd -P)"
CURTIME=$(date)
echo "$CURTIME: starting HyPhy BS-REL." 

#do run
cd tree1
if [ -s $TR1IN ];
then
	HYPHYMP BASEPATH=/n/home12/tsackton/sw/source/hyphy/res/TemplateBatchFiles/ CPU=1 BranchSiteREL.bf < input > output
fi

CURTIME=$(date)
echo "$CURTIME: Finished Tree1."

cd ..
cd tree2
if [ -s $TR2IN ]
then
	HYPHYMP BASEPATH=/n/home12/tsackton/sw/source/hyphy/res/TemplateBatchFiles/ CPU=1 BranchSiteREL.bf < input > output
fi

CURTIME=$(date)
echo "$CURTIME: Finished Tree2."

cd ..
