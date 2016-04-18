#post process moa results to prepare for R analysis

for NEUT in neut_ver1 neut_ver2 neut_ver3
do
#	for SAMP in allRatite.out anoDid.out aptHaa.out aptOwe.out aptRow.out casCas.out Casuar.out droNov.out Kiwi.out Rhea.out rheAme.out rhePen.out strCam.out Tinamou.out
	for SAMP in basalPaleo.out
	do
		
		head -n 1 batch1/$NEUT/$SAMP | cut -f4-9 > header.temp
		tail -n +2 batch1/$NEUT/$SAMP > ${SAMP}_${NEUT}.1.temp
		tail -n +2 batch2/$NEUT/$SAMP > ${SAMP}_${NEUT}.2.temp
		tail -n +2 batch3/$NEUT/$SAMP > ${SAMP}_${NEUT}.3.temp
		tail -n +2 batch4/$NEUT/$SAMP > ${SAMP}_${NEUT}.4.temp
		cat ${SAMP}_${NEUT}.1.temp ${SAMP}_${NEUT}.2.temp ${SAMP}_${NEUT}.3.temp ${SAMP}_${NEUT}.4.temp | cut -f4-9 | awk '!seen[$0]++' > ${SAMP}_${NEUT}.temp
		cat header.temp ${SAMP}_${NEUT}.temp > ${SAMP}_${NEUT}.results
		rm *.temp
	done
done
