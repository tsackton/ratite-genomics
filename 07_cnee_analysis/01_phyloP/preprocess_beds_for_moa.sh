#code to preprocess beds in order to add moa sequence to alignments for phyloP processing
#start with final_cnees_long.bed which is just the long (>= 50bp) CNEEs from the initial processing

#liftover to all species
#working in /n/regal/edwards_lab/ratites/wga/phast/moa/cnee_liftover
for SP in allMis aptFor balReg chaVoc corBra droNov lepDis nipNip pygAde taeGut allSin aptHaa calAnn cheMyd croPor eudEle fulGla melGal notPer rheAme tinGut anaPla aptOwe casCas chrPic cryCin falPer gavGan melUnd picPub rhePen anoCar aptRow chaPel colLiv cucCan ficAlb halLeu mesUni pseHum strCam
do
	halLiftover --outPSLWithName ~/ratite_scratch/wga/ratite_final_20150627/ratiteAlign.hal galGal final_cnees_long.bed $SP $SP.psl &
done

#now get 4d sites, working in /n/regal/edwards_lab/ratites/wga/phast/moa/neutMods_withMoa
awk 'NF==12{print}{}' galGal4_4d.bed > galGal4_4d_bed12.bed
bedtools bed12tobed6 -i galGal4_4d_bed12.bed > galGal4_4d_bed6.bed 

#replace ids
awk 'BEGIN {FS="\t"; OFS="\t"} {$4="fourD."NR; print}' galGal4_4d_bed6.bed > galGal4_4d_bed6_withid.bed 

#now lift over the bed to each species
for SP in allMis aptFor balReg chaVoc corBra droNov lepDis nipNip pygAde taeGut allSin aptHaa calAnn cheMyd croPor eudEle fulGla melGal notPer rheAme tinGut anaPla aptOwe casCas chrPic cryCin falPer gavGan melUnd picPub rhePen anoCar aptRow chaPel colLiv cucCan ficAlb halLeu mesUni pseHum strCam
do
	halLiftover ~/ratite_scratch/wga/ratite_final_20150627/ratiteAlign.hal galGal galGal4_4d_bed6_withid.bed $SP $SP.4d.bed &
done
