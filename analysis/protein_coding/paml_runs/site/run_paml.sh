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


#clade models 
for TREETYPE in spt gt
do
	if [[ $TREETYPE == "spt" ]]
	then
		if [[ $MAXSP -gt 1 ]]
		then
			continue;
		else		
			#named clades
			RHEA=$(comm -12 <(echo -e "rheAme\nrhePen\nrheAme-rhePen\nrhePen-rheAme" | sort) <(nw_labels $HOG.final_$TREETYPE.nwk | sort))
			KIWI=$(comm -12 <(echo -e "aptHaa\naptRow\naptOwe\naptHaa-aptRow\naptHaa-aptOwe\naptRow-aptOwe\naptRow-aptHaa\naptOwe-aptRow\naptOwe-aptHaa" | sort) <(nw_labels $HOG.final_$TREETYPE.nwk | sort))
			ECAS=$(comm -12 <(echo -e "droNov\ncasCas\ndroNov-casCas\ncasCas-droNov\n" | sort) <(nw_labels $HOG.final_$TREETYPE.nwk | sort))
			OST=$(comm -12 <(echo -e "strCam") <(nw_labels $HOG.final_$TREETYPE.nwk | sort))
		fi
	else
		if [[ $MAXSP -gt 1 ]]
		then
			RHEA=$(grep 'rheAme\|rhePen' $HOG.sp.key | cut -f1,1)
			ECAS=$(grep 'droNov\|casCas' $HOG.sp.key | cut -f1,1)
			KIWI=$(grep 'aptHaa\|aptRow\|aptOwe' $HOG.sp.key | cut -f1,1)
			OST=$(grep 'strCam' $HOG.sp.key | cut -f1,1)
		else
                        RHEA=$(grep 'rheAme\|rhePen' $HOG.sp.key | cut -f2,2)
                        ECAS=$(grep 'droNov\|casCas' $HOG.sp.key | cut -f2,2)
                        KIWI=$(grep 'aptHaa\|aptRow\|aptOwe' $HOG.sp.key | cut -f2,2)
			OST=$(grep 'strCam' $HOG.sp.key | cut -f2,2)
		fi	
	fi
	
	RATITE="$RHEA $KIWI $ECAS $OST"	
	
	#label branches with tree_doctor
	tree_doctor -n -N -l $(echo $RHEA | awk -v OFS="," '$1=$1'):1 $HOG.final_$TREETYPE.nwk > rhea.temp
	tree_doctor -n -N -l $(echo $RATITE | awk -v OFS="," '$1=$1'):1 $HOG.final_$TREETYPE.nwk > ratite.temp
	tree_doctor -n -N -l $(echo $KIWI | awk -v OFS="," '$1=$1'):1 $HOG.final_$TREETYPE.nwk > kiwi.temp
	tree_doctor -n -N -l $(echo $ECAS | awk -v OFS="," '$1=$1'):1 $HOG.final_$TREETYPE.nwk > ecas.temp
	tree_doctor -n -N -l $(echo $OST | awk -v OFS="," '$1=$1'):1 $HOG.final_$TREETYPE.nwk > ost.temp

	cat *.temp >> $HOG.final_clade_$TREETYPE.nwk;
	rm *.temp
	
done

#concatenate final trees
if [[ $MAXSP -eq 1 ]]
then
	perl -p -i -e 's/\w+-\w+//g' $HOG.final_spt.nwk
	perl -p	-i -e 's/\w+-\w+//g' $HOG.final_clade_spt.nwk
	cat $HOG.final_spt.nwk $HOG.final_gt.nwk | wc -l > $HOG.final.nwk
	cat $HOG.final_clade_spt.nwk $HOG.final_clade_gt.nwk | wc -l > $HOG.finalclade.nwk
	cat $HOG.final_spt.nwk $HOG.final_gt.nwk >> $HOG.final.nwk
	cat $HOG.final_clade_spt.nwk $HOG.final_clade_gt.nwk >> $HOG.finalclade.nwk
else
	cat $HOG.final_gt.nwk | wc -l > $HOG.final.nwk
	cat $HOG.final_clade_gt.nwk | wc -l > $HOG.finalclade.nwk
	cat $HOG.final_gt.nwk >> $HOG.final.nwk
	cat $HOG.final_clade_gt.nwk >> $HOG.finalclade.nwk
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

for BASECTL in $(ls $WORKDIR/codeml.*.ctl)
do
	CTL=${BASECTL##*/}
	cp $BASECTL $HOG.$CTL
	MODEL=${CTL##codeml.}
	sed -i "s/SEQINPUT/$HOG.phy/" $HOG.$CTL
	if [[ $CTL == *"branch"* ]] 
	then sed -i "s/TREEINPUT/$HOG.finalclade.nwk/" $HOG.$CTL;
	else sed -i "s/TREEINPUT/$HOG.final.nwk/" $HOG.$CTL;
	fi
	sed -i "s/OUTPUT/${MODEL%%.ctl}.out/" $HOG.$CTL
done

CURTIME=$(date)
echo "$CURTIME: Starting PAML runs."

for RUN in $(ls *.ctl)
do
	mkdir -p $RUN.out
	if [[ -e $RUN.out/DONE ]]
	then
		echo "NOTE: Skipping completed run $RUN"
		continue
	fi
	if [[ -e $HOG.codeml.branch.ctl.out/DONE && $RUN == "$HOG.codeml.br.ctl" ]]
	then
		echo "NOTE: skipping completed branch run for $RUN"
	fi 
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
