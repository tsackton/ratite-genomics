#!/bin/bash

##PAML SCRIPT / PIPELINE##

#get info, write to log

WORKDIR=$(pwd -P)
ALNDIR=/n/holylfs/LABS/edwards_lab/Users/tsackton/ratites/homology/oma_nucl_aligns
TREEDIR=/n/holylfs/LABS/edwards_lab/Users/tsackton/ratites/homology/oma_raxml_trees
HOG=$1
CURTIME=$(date)

echo "$CURTIME: Running PAML for $HOG"
INPUT="$ALNDIR/HOG2_$HOG.fa"
GENETREE="$TREEDIR/RAxML_bipartitions.HOG2_$HOG"
SEQNUM=$(grep -c ">" $INPUT)
#GAPTHRESH=$(bc -l <<< "1-(2/$SEQNUM)")
GAPTHRESH=0.80

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
nw_topology ${HOG}_nbs.orig.nwk > $HOG.init_gt.nwk

#copy input
cp $INPUT $HOG.init.fa
perl -p -i -e 's/[=]//g' $HOG.init.fa
perl -p -i -e 's/[=]//g' $HOG.init_gt.nwk

#trim initial fasta
trimal -in $HOG.init.fa -out $HOG.trimmed.fa -fasta -gt $GAPTHRESH

#make species key
grep "^>" $HOG.trimmed.fa | perl -p -e 's/>([A-Za-z]*)(_.*)/$1$2\t$1/' > $HOG.sp.key

#prune gene tree
GENELIST=$(cut -f1,1 $HOG.sp.key)
nw_prune -v $HOG.init_gt.nwk $GENELIST > $HOG.final_gt.nwk

#make tree files -- complicated bit
if [[ $MAXSP -gt 1 ]] 
then
	cp $HOG.trimmed.fa $HOG.fa
else
	#this means we first rename species 
	perl -p -e 's/>([A-Za-z]*)_.*/>$1/' $HOG.trimmed.fa > $HOG.fa
	nw_rename $HOG.final_gt.nwk $HOG.sp.key > $HOG.final_gt.renamed.nwk
	mv $HOG.final_gt.renamed.nwk $HOG.final_gt.nwk
	SPLIST=$(cut -f2,2 $HOG.sp.key)
	nw_prune -v $WORKDIR/oma_tree.nh $SPLIST > $HOG.final_spt.nwk
	tree_doctor -n -a $HOG.final_spt.nwk > $HOG.final_spt_named.nwk
	mv $HOG.final_spt.nwk $HOG.spt_nonames.nwk
	mv $HOG.final_spt_named.nwk $HOG.final_spt.nwk	
fi

#concatenate final trees
if [[ $MAXSP -eq 1 ]]
then
	perl -p -i -e 's/\w+-\w+//g' $HOG.final_spt.nwk
	cat $HOG.final_spt.nwk $HOG.final_gt.nwk | wc -l > $HOG.final.nwk
	cat $HOG.final_spt.nwk $HOG.final_gt.nwk >> $HOG.final.nwk
else
	cat $HOG.final_gt.nwk | wc -l > $HOG.final.nwk
	cat $HOG.final_gt.nwk >> $HOG.final.nwk
fi

#trim final fasta
trimal -in $HOG.fa -out $HOG.phy -phylip_paml -gt $GAPTHRESH

#all inputs present
echo "All inputs created."
CURTIME=$(date)
echo "$CURTIME: Starting PAML setup."

#now just need to sort out running options, based on codeml files
#codeml.ancrec.ctl --> requires bifurcating species tree, run only for non-dups
#codeml.branch.ctl --> free ratio dN/dS, run for all alignments
#codeml.branchsite.ctl --> testing positive selection on ratite branches, need to run many trees
#codeml.branchclade.ctl --> testing different rates on different clades. NOT RUN.
#codeml.site.ctl --> basic site test

#set up paml runs -- TO DO##

for BASECTL in $(ls $WORKDIR/codeml.ancrec.ctl)
do
	CTL=${BASECTL##*/}
	cp $BASECTL $HOG.$CTL
	MODEL=${CTL##codeml.}
	sed -i "s/SEQINPUT/$HOG.phy/" $HOG.$CTL
	sed -i "s/OUTPUT/${MODEL%%.ctl}.out/" $HOG.$CTL
	sed -i "s/TREEINPUT/$HOG.final.nwk/" $HOG.$CTL
done

CURTIME=$(date)
echo "$CURTIME: Starting PAML runs."

for RUN in $(ls *.ctl)
do
	mkdir -p $RUN.out
	mv $RUN $RUN.out/
	cp $HOG.phy $RUN.out
	cp $HOG*.nwk $RUN.out
	cd $RUN.out
	codeml $RUN
	CURTIME=$(date)
	echo "$CURTIME: Finished $RUN"
	touch DONE
	cd ..
done
