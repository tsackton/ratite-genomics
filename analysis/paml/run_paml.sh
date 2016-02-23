#!/bin/bash

##PAML SCRIPT / PIPELINE##

#first trim with trimal, remove everything except 

WORKDIR=$(pwd -P)
ALNDIR=/n/regal/edwards_lab/ratites/align_oma/final_prank/alignments/prank/ogs_aligned/dna
HOG=$1
CURTIME=$(date)

echo "$CURTIME: Running PAML for $HOG"
INPUT="$ALNDIR/$HOG.dna.best.fas"
SEQNUM=$(grep -c ">" $INPUT)
GAPTHRESH=$(bc -l <<< "1-(2/$SEQNUM)")

#setup working dir
SUBDIRNUM=`expr ${HOG##HOG} % 100`
printf -v SUBDIR "%04d" $SUBDIRNUM
PAMLDIR=$SUBDIR/$HOG	
mkdir -p $PAMLDIR

cd $PAMLDIR

#trim
trimal -in $INPUT -out $HOG.phy -phylip_paml -gt $GAPTHRESH

#now prune tree
SPLIST=$(grep ">" $INPUT | sed 's/>//' | paste -sd ' ' -)
nw_prune -v $WORKDIR/prank_tree.nwk $SPLIST > $HOG.nwk

#now make branch model tree with two rates (ratites vs others)
#just use flightless tips - this is conservative as some fraction of internal nodes likely represent flightless phenotypes
cp $HOG.nwk $HOG.clade.nwk
for SPTOFIX in droNov casCas aptHaa strCam rheAme rhePen
do
	sed -i "s/$SPTOFIX/$SPTOFIX #1/" $HOG.clade.nwk
done

#codeml files:
#codeml.site.ctl <- NOT RUN
#codeml.ancrec.ctl <- ancestral reconstruction with M0
#codeml.branch.ctl <- branch model
#codeml.branchsite.ctl <- branch/site model
#codeml.branchclade.ctl <- clade model

for BASECTL in $(ls $WORKDIR/codeml.*.ctl)
do
	CTL=${BASECTL##*/}
	cp $BASECTL $HOG.$CTL
	MODEL=${CTL##codeml.}
	sed -i "s/SEQINPUT/$HOG.phy/" $HOG.$CTL
	if [[ $CTL == *"branch"* ]] 
	then sed -i "s/TREEINPUT/$HOG.clade.nwk/" $HOG.$CTL;
	else sed -i "s/TREEINPUT/$HOG.nwk/" $HOG.$CTL;
	fi
	sed -i "s/OUTPUT/${MODEL%%.ctl}.out/" $HOG.$CTL
done

for RUN in $(ls *.ctl)
do
	mkdir $RUN.out
	mv $RUN $RUN.out/
	cp $HOG.phy $RUN.out
	cp $HOG*.nwk $RUN.out
	cd $RUN.out
	codeml $RUN
	CURTIME=$(date)
	echo "$CURTIME: Finished $RUN"
	cd ..
done
