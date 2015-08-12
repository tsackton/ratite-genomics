#!/bin/bash

#process junctions.bed to keep only high-quality junctions (score > 5) and then convert to gffs
#merge alt-species into one GFF
#keep same-species (droNov/droNov and notPer/notPer) independent

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

#kiwi, cassowary, rheas have just altest of both kiwi and notper
gff3_merge -o aptHaa.altest.gff aptHaa*.junctions.bed.gff
gff3_merge -o aptRow.altest.gff aptRow*.junctions.bed.gff
gff3_merge -o aptOwe.altest.gff aptOwe*.junctions.bed.gff
gff3_merge -o casCas.altest.gff casCas*.junctions.bed.gff
gff3_merge -o rheAme.altest.gff rheAme*.junctions.bed.gff
gff3_merge -o rhePen.altest.gff rhePen*.junctions.bed.gff

