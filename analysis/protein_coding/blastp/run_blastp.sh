module load blast/2.2.29+-fasrc01
for SP in *.phr
do
	TARGET=${SP%%.*}
	blastp -query Gallus_gallus.Galgal4.pep.all.fa -db $TARGET -out ${TARGET}_galGal_blastp.out -outfmt 6 -evalue 1e-3 -num_threads 16
done
