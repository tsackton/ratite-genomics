module load blast/2.2.29+-fasrc01
for FASTA in *.fa*;
do
	makeblastdb -in $FASTA -dbtype 'prot' -out ${FASTA%%.*}
done

