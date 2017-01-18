for BLASTP in *_blastp.out
do
	sort -u -k1,1 $BLASTP | awk -v name="$BLASTP" 'BEGIN{OFS="\t"} {print name, $0}' > $BLASTP.besthit &
done
