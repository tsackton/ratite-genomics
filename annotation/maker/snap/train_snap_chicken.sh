#!/bin/bash

#code to do initial SNAP training on chicken

#1. train SNAP on chicken: (based on https://biowize.wordpress.com/2012/06/01/training-the-snap-ab-initio-gene-predictor/ and SNAP readme)
#a. Get chicken genome, chromosomes only:

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 Z
do
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000002315.3_Gallus_gallus-4.0/GCF_000002315.3_Gallus_gallus-4.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr$CHR.fna.gz
done
cat *.gz > galgal4chr.fa.gz
gunzip galgal4chr.fa.gz
rm *.gz
	
#b. Get chicken GFF from NCBI
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000002315.3_Gallus_gallus-4.0/GCF_000002315.3_Gallus_gallus-4.0_genomic.gff.gz
gunzip GCF_000002315.3_Gallus_gallus-4.0_genomic.gff.gz
	
#c. Filter chicken GFF 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000002315.3_Gallus_gallus-4.0/GCF_000002315.3_Gallus_gallus-4.0_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc
for ACC in $(egrep "^[0-9Z]" chr2acc | grep -v "^32" | cut -f2,2)
do
	egrep "^$ACC" GCF_000002315.3_Gallus_gallus-4.0_genomic.gff >> galgal4_chr.gff
done
grep "BestRefSeq" galgal4_chr.gff > galgal4_filt_chr.gff
grep -v "exception=" galgal4_filt_chr.gff > galgal4_filt_chr_good.gff

#use gffread to remove cds with stop codons and other problems
module load cufflinks
gffread galgal4_filt_chr_good.gff -g galgal4chr.fa -o galgal4_final.gff -V -C 

#merge fasta info to end of gff
echo "##FASTA" >> galgal4_final.gff
cat galgal4chr.fa >> galgal4_final.gff	

#d. convert gff to zff
maker2zff -n galgal4_final.gff
		
#e. check quality of zff files
fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
#MODEL589 skipped due to errors
#MODEL1127 1 1 4 - errors(2): gene:misordered_Einit gene:misordered_Eterm
#MODEL1127 skipped due to errors
#MODEL2463 1 1 12 - errors(2): gene:misordered_Einit gene:misordered_Eterm
#MODEL2463 skipped due to errors

fathom genome.ann genome.dna -validate > validate.log 2>&1
grep -v "MODEL589" genome.ann | grep -v "MODEL1127" | grep -v "MODEL2463" > genome.fixed.ann
mv genome.ann genome.ann.old
mv genome.fixed.ann genome.ann
fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
	
#f. run SNAP training algo
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
hmm-assembler.pl chicken params/ > chicken.hmm
