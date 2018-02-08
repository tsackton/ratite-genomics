#!/bin/bash

##PAML SCRIPT / PIPELINE##

#get info, write to log
#move to correct dir

WORKDIR=$(pwd -P)
HOG=$1
SUBDIRNUM=`expr ${HOG##HOG} % 100`
printf -v SUBDIR "%04d" $SUBDIRNUM
PAMLDIR=$SUBDIR/$HOG	
cd $PAMLDIR

#set up paml runs

for BASECTL in $(ls $WORKDIR/codeml.br.ctl)
do
	CTL=${BASECTL##*/}
	cp $BASECTL $HOG.$CTL
	MODEL=${CTL##codeml.}
	sed -i "s/SEQINPUT/$HOG.phy/" $HOG.$CTL
	sed -i "s/TREEINPUT/$HOG.final.nwk/" $HOG.$CTL;
	sed -i "s/OUTPUT/${MODEL%%.ctl}.out/" $HOG.$CTL
done

CURTIME=$(date)
echo "$CURTIME: Starting PAML runs."

for RUN in $(ls *.br.ctl)
do
	mkdir -p $RUN.out
	cp $RUN $RUN.out
	cp $HOG.phy $RUN.out
	cp $HOG*.nwk $RUN.out
	cd $RUN.out
	echo "$(pwd -P)"
	codeml $RUN
	CURTIME=$(date)
	echo "$CURTIME: Finished $RUN"
	touch DONE
	cd ..
done
