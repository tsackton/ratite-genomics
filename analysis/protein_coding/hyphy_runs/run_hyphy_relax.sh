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
mkdir -p tree1/all
mkdir -p tree1/tips
mkdir -p tree2/all
mkdir -p tree2/tips
ALIGN=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/DATA/$RUNDIR/$HOG.phy
INTREE1=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/DATA/$RUNDIR/$HOG.tree1.nwk
INTREE2=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/DATA/$RUNDIR/$HOG.tree2.nwk
OUTTREE1=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/RELAX_Ratite/$RUNDIR/tree1/$HOG.relax.nwk
OUTTREE2=/n/holylfs/LABS/edwards_lab/tsackton/HyPhy/RELAX_Ratite/$RUNDIR/tree2/$HOG.relax.nwk
python $WORKDIR/make_relax_trees.py $INTREE1 $OUTTREE1
python $WORKDIR/make_relax_trees.py $INTREE2 $OUTTREE2
#UGLY!
perl -p -i -e 's/\[&&NHX:mark=//g' $OUTTREE1 
perl -p -i -e 's/\]//g' $OUTTREE1
perl -p -i -e 's/\[&&NHX:mark=//g' $OUTTREE2            
perl -p -i -e 's/\]//g' $OUTTREE2 
perl -p -e 's/(Leaf)|(Internal)//g' $OUTTREE1 > $OUTTREE1.all.nwk
perl -p -e 's/(Leaf)|({RatiteInternal})//g' $OUTTREE1 > $OUTTREE1.tips.nwk
perl -p -e 's/(Leaf)|(Internal)//g' $OUTTREE2 > $OUTTREE2.all.nwk       
perl -p -e 's/(Leaf)|({RatiteInternal})//g' $OUTTREE2 > $OUTTREE2.tips.nwk  
echo -e "1\n$ALIGN\n$OUTTREE1.tips.nwk\n1\n2\n" > tree1/tips/input
echo -e "1\n$ALIGN\n$OUTTREE1.all.nwk\n1\n2\n" > tree1/all/input
echo -e "1\n$ALIGN\n$OUTTREE2.tips.nwk\n1\n2\n" > tree2/tips/input
echo -e "1\n$ALIGN\n$OUTTREE2.all.nwk\n1\n2\n" > tree2/all/input

echo "$(pwd -P)"
CURTIME=$(date)
echo "$CURTIME: starting HyPhy RELAX." 
for TREEVER in tree1 tree2
    do
    for ATYPE in all tips
         do
         cd $TREEVER/$ATYPE
         echo "Running...."
         HYPHYMP BASEPATH=/n/home12/tsackton/sw/source/hyphy/res/TemplateBatchFiles/ CPU=1 RELAX.bf < input > output
         echo "$CURTIME: Finished $TREEVER, $ATYPE" 
         cd ../..
    done
done
