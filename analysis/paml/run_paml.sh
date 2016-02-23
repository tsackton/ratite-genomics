#!/bin/bash

##PAML SCRIPT / PIPELINE##

#get info, write to log

WORKDIR=$(pwd -P)
ALNDIR=/n/edwards_lab/Users/tsackton/ratites/homology/oma_nucl_aligns
TREEDIR=/n/edwards_lab/Users/tsackton/ratites/homology/oma_raxml_trees
HOG=$1
CURTIME=$(date)

echo "$CURTIME: Running PAML for $HOG"
INPUT="$ALNDIR/HOG2_$HOG.fa"
GENETREE="$TREEDIR/RAxML_bipartitions.HOG2_$HOG"
SEQNUM=$(grep -c ">" $INPUT)
GAPTHRESH=$(bc -l <<< "1-(2/$SEQNUM)")

#setup working dir
SUBDIRNUM=`expr ${HOG##HOG} % 100`
printf -v SUBDIR "%04d" $SUBDIRNUM
PAMLDIR=$SUBDIR/$HOG	
mkdir -p $PAMLDIR

cd $PAMLDIR

#detect if duplicates are present
MAXSP=$(grep ">" $INPUT| perl -p -e 's/>([A-Za-z]*)_.*/$1/' | sort | uniq -c | awk -F" " '{print $1}' | sort -nr | head -n 1)

echo "Found max of $MAXSP sequences per species."

#now set up tree -- make two versions, one that is based on oma tree, one that is based on gene tree
#first copy gene tree and process with TreeCollapseCL4
cp $GENETREE $HOG.orig.nwk
java -jar /n/home12/tsackton/sw/bin/TreeCollapseCL4.jar -f $HOG.orig.nwk -b 75 -t O -nbs
nw_topology ${HOG}_75coll.orig.nwk > $HOG.final_gt.nwk

#make species key
grep "^>" $INPUT | perl -p -e 's/>([A-Za-z]*)(_.*)/$1$2\t$1/' > $HOG.sp.key

#first need to rename sequences in alignment if there are no duplications
if [[ $MAXSP > 1 ]] 
then
	#this means we have duplications and so will need to just use the gene tree
	cp $INPUT $HOG.fa
else
	#this means we first rename species 
	perl -p -e 's/>([A-Za-z]*)_.*/>$1/' $INPUT > $HOG.fa
	nw_rename $HOG.final_gt.nwk $HOG.sp.key > $HOG.final_gt.renamed.nwk
	mv $HOG.final_gt.renamed.nwk $HOG.final_gt.nwk
	SPLIST=$(cut -f2,2 $HOG.sp.key)
	nw_prune -v $WORKDIR/oma_tree.nh $SPLIST > $HOG.final_spt.nwk
fi

#have two tree files for single-copy genes, one for duplicated genes
#make phylip
#trim
trimal -in $HOG.fa -out $HOG.phy -phylip_paml -gt $GAPTHRESH

#all inputs present
echo "All inputs created."
CURTIME=$(date)
echo "$CURTIME: Starting PAML setup."

#now just need to sort out running options, based on codeml files
#codeml.ancrec.ctl --> requires bifurcating species tree, run only for non-dups
#codeml.branch.ctl --> free ratio dN/dS, run for all alignments
#codeml.branchsite.ctl --> testing positive selection on ratite branches, need to run many trees
#codeml.branchclade.ctl --> testing different rates on different clades
#codeml.site.ctl --> basic site test


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
