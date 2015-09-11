#!/bin/bash

#process junctions.bed to keep only high-quality junctions (score > 5) and then convert to gffs
#merge alt-species into one GFF
#keep same-species (droNov/droNov and notPer/notPer) independent

mkdir output_for_maker
cd output_for_maker

#first filter with awk
for BED in $(ls ../*/junctions.bed);
do
	DIR=${BED%/*}
	TARGET=${DIR#*/}
	awk '{if($5 >= 5) print $0}' $BED > $TARGET.junctions.bed
done

#convert to GFF
for BED in 	$(ls *.bed)
do
	tophat2gff3 $BED > $BED.gff
done

#merge and rename

#dronov has both
cp droNov.dronov.junctions.bed.gff droNov.est.gff
cp droNov.kiwi.junctions.bed.gff droNov.altest.gff

#notper has just est
cp notPer.notper.junctions.bed.gff notPer.est.gff

#tinamous just have altest
cp cryCin.notper.junctions.bed.gff cryCin.altest.gff
cp eudEle.notper.junctions.bed.gff eudEle.altest.gff

#kiwi have both
cp aptHaa.kiwi.junctions.bed.gff aptHaa.est.gff
cp aptHaa.dronov.junctions.bed.gff aptHaa.altest.gff
cp aptRow.kiwi.junctions.bed.gff aptRow.est.gff
cp aptRow.dronov.junctions.bed.gff aptRow.altest.gff
cp aptOwe.kiwi.junctions.bed.gff aptOwe.est.gff
cp aptOwe.dronov.junctions.bed.gff aptOwe.altest.gff

#cassowary, rheas have just altest of both kiwi and notper
#added -l option to deal with multiple indentical ids
gff3_merge -o casCas.altest.gff casCas*.junctions.bed.gff
gff3_merge -l -o rheAme.altest.gff rheAme*.junctions.bed.gff
gff3_merge -l -o rhePen.altest.gff rhePen*.junctions.bed.gff

