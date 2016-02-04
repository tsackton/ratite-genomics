##NOTE: code below uses bedtools v2.24.0

#define a bash division function
divide() {
   echo "scale=25;$1/$2" | bc
}

#get lowe et al CNEEs in galGal4, NCBI style coordinates
#step one is liftover, then replace chrs using replace chrs script; info in galGalAnnot directory
cp ~/ratite_store/ratite-genomics/annotation/galGalAnnot/LoweCNEE/lowe_cnees.bed .
cp ../phastCons/most_conserved_*.bed .
cp ../replace_chrs.pl .

#get galGal exons using Phil's GFF_CDS script and GFF_BED script, then sort and merge:
python ~/ratite_store/ratite-genomics/analysis/missing_genes/rna_intersect/GFF_CDS_parser.py ../GCF_000002315.3_Gallus_gallus-4.0_genomic.gff -o galGal_exons.gff -t exon
python ~/ratite_store/ratite-genomics/analysis/missing_genes/rna_intersect/GFF_CDS_parser.py ../GCF_000002315.3_Gallus_gallus-4.0_genomic.gff -o galGal_cds.gff -t CDS
python ~/ratite_store/ratite-genomics/analysis/missing_genes/rna_intersect/GFF_CDS_parser.py ../GCF_000002315.3_Gallus_gallus-4.0_genomic.gff -o galGal_gene.gff -t gene
python ~/ratite_store/ratite-genomics/analysis/missing_genes/rna_intersect/Convert_GFF_to_BED.py galGal_exons.gff -o galGal_exons.bed
python ~/ratite_store/ratite-genomics/analysis/missing_genes/rna_intersect/Convert_GFF_to_BED.py galGal_cds.gff -o galGal_CDS.bed
python ~/ratite_store/ratite-genomics/analysis/missing_genes/rna_intersect/Convert_GFF_to_BED.py galGal_gene.gff -o galGal_gene.bed

sort -k1,1 -k2,2n -k3,3n -u galGal_exons.bed | bedtools merge -i - -s -d -1 -c 4 -o distinct > chicken_exons.bed
sort -k1,1 -k2,2n -k3,3n -u galGal_CDS.bed | bedtools merge -i - -s -d -1 -c 4 -o distinct > chicken_CDS.bed
bedtools subtract -s -a chicken_exons.bed -b chicken_CDS.bed > chicken_nonCDS.bed

#repeat for ensembl genes
#download gff
wget ftp://ftp.ensembl.org/pub/release-83/gff3/gallus_gallus/Gallus_gallus.Galgal4.83.gff3.gz
gunzip Gallus_gallus.Galgal4.83.gff3.gz
#download assembly report
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000002315.3_Gallus_gallus-4.0/GCF_000002315.3_Gallus_gallus-4.0_assembly_report.txt
#make conversion file
grep -v "^#" GCF_000002315.3_Gallus_gallus-4.0_assembly_report.txt | grep "scaffold" | cut -f5,7 > ensToNCBI.chrlist
grep -v "^#" GCF_000002315.3_Gallus_gallus-4.0_assembly_report.txt | grep "assembled" | cut -f3,7 >> ensToNCBI.chrlist 
./replace_chrs.pl Gallus_gallus.Galgal4.83.gff3 
#get genes
python ~/ratite_store/ratite-genomics/analysis/missing_genes/rna_intersect/GFF_CDS_parser.py  Gallus_gallus.Galgal4.83.gff3.fixed -o galGal_gene_ens.gff -t gene
python ~/ratite_store/ratite-genomics/analysis/missing_genes/rna_intersect/Convert_GFF_to_BED.py galGal_gene_ens.gff -o galGal_gene_ens.bed

#check for coverage in CDS, non-CDS, cnees, genome
bedtools coverage -a chicken_CDS.bed -b most_conserved_final.tree2.bed > CDS.coverage.tree2
awk '{olen += $6; sum+=$7}END{print olen/sum}' CDS.coverage.tree2
bedtools coverage -a chicken_nonCDS.bed -b most_conserved_final.tree2.bed > nonCDS.coverage.tree2
awk '{olen += $6; sum+=$7}END{print olen/sum}' nonCDS.coverage.tree2
bedtools coverage -a lowe_cnees.bed -b most_conserved_final.tree2.bed > cnee.coverage.tree2
awk '{olen += $6; sum+=$7}END{print olen/sum}' cnee.coverage.tree2
#genome cov
GSIZE=$(awk '{sum+=$2}END{print sum}' ../galGal.chromsizes)
CLEN=$(awk '{sum+=$3-$2}END{print sum}' most_conserved_final.tree2.bed)
divide $CLEN $GSIZE

#merge nearby ces in existing bed files
bedtools merge -i most_conserved_final.tree2.bed -d 5 -c 4 -o collapse -delim "|" > merged_ces.tree2.bed

#recheck coverage
bedtools coverage -a chicken_CDS.bed -b merged_ces.tree2.bed > CDS.coverage.merged.tree2
awk '{olen += $6; sum+=$7}END{print olen/sum}' CDS.coverage.merged.tree2
CLEN=$(awk '{sum+=$3-$2}END{print sum}' merged_ces.tree2.bed)
divide $CLEN $GSIZE

#get counts of elements per cds
bedtools intersect -a chicken_CDS.bed -b most_conserved_final.tree2.bed  -c > tree2.CDS.count
bedtools intersect -a chicken_CDS.bed -b merged_ces.tree2.bed  -c > tree2.merged.CDS.count

#analyze in R
merge<-read.table("tree2.merged.CDS.count", header=F, sep="\t")
init<-read.table("tree2.CDS.count", header=F, sep="\t")
sum(init$V5==1)/length(init$V5)
#[1] 0.2183828
sum(merge$V5==1)/length(merge$V5)
#[1] 0.3813054

wc -l merged_ces.tree2.bed 
#1956236 merged_ces.tree2.bed
wc -l most_conserved_final.tree2.bed 
#2299364 most_conserved_final.tree2.bed

#5 bp is arbitrary but looks pretty good so we'll go with that
bedtools merge -i most_conserved_final.tree1.bed -d 5 -c 4 -o collapse -delim "|" > merged_ces.tree1.bed
bedtools merge -i most_conserved_noratite_final.tree2.bed -d 5 -c 4 -o collapse -delim "|" > merged_ces.nr.tree2.bed

#get final overlaps
for BED in merged_ces* 
do
	bedtools intersect -a chicken_CDS.bed -b $BED -c > $BED.CDS.count
	bedtools intersect -a lowe_cnees.bed -b $BED -c > $BED.lowe_cnee.count
	bedtools intersect -a chicken_nonCDS.bed -b $BED -c > $BED.nonCDS.count
done

#rename / renumber
awk 'BEGIN {FS="\t"; OFS="\t"} {$4="mCE"NR; print}' merged_ces.tree1.bed > final_ces.tree1.bed
awk 'BEGIN {FS="\t"; OFS="\t"} {$4="mCE"NR; print}' merged_ces.tree2.bed > final_ces.tree2.bed
awk 'BEGIN {FS="\t"; OFS="\t"} {$4="mnrCE"NR; print}' merged_ces.nr.tree2.bed > final_ces_noratite.tree2.bed

#make a "intersection" file that is just the consistent conserved elements between tree1 and tree2 for phylogenetic analysis:
bedtools intersect -a final_ces.tree2.bed -b final_ces.tree1.bed -f 0.50 -r -u > final_ces.intersection.bed

#sort gene files
sort -k1,1 -k2,2n galGal_gene.bed > galGal_ncbi.sorted.bed 
perl -p -i -e 's/,Start:\d+,Stop:\d+,Strand:.//' galGal_ncbi.sorted.bed
perl -p -i -e 's/CGNC:\d+,//' galGal_ncbi.sorted.bed
perl -p -i -e 's/GeneID://' galGal_ncbi.sorted.bed
sort -k1,1 -k2,2n galGal_gene_ens.bed > galGal_ens.sorted.bed 
perl -p -i -e 's/ /_/g' galGal_ens.sorted.bed 
perl -p -i -e 's/ID=gene:(\w+)\S+/$1/' galGal_ens.sorted.bed

for BED in final_ces*.tree?.bed
do
	#get counts for each feature
	bedtools intersect -b $BED -a chicken_CDS.bed -c > $BED.ct.CDS
	bedtools intersect -b $BED -a chicken_exons.bed -c > $BED.ct.exon
	bedtools intersect -b $BED -a galGal_gene.bed -c > $BED.ct.gene
	bedtools intersect -a lowe_cnees.bed -b $BED -c > $BED.ct.lowe	
	#get overlaps
	bedtools intersect -a $BED -b chicken_CDS.bed -c > $BED.annot.CDS
	bedtools intersect -a $BED -b chicken_exons.bed -c > $BED.annot.exon
	bedtools intersect -a $BED -b galGal_gene.bed -c > $BED.annot.gene
	#get lowe cnee overlap id
	bedtools intersect -a $BED -b lowe_cnees.bed -loj > $BED.annot.lowe
	#get closest gene - NCBI
	bedtools closest -a $BED -b galGal_ncbi.sorted.bed  -D "b" -t "all" > $BED.annot.closest_genes_ncbi
	#get closest gene - Ensembl
	bedtools closest -a $BED -b galGal_ens.sorted.bed -D "b" -t "all" > $BED.annot.closest_genes_ens
done

#get closest genes defined by approximate TSS
./get_approx_TSS.sh galGal_ens.sorted.bed > galGal_ens_approxTSS.bed
./get_approx_TSS.sh galGal_ncbi.sorted.bed > galGal_ncbi_approxTSS.bed

for BED in final_ces*.tree?.bed
do
	bedtools closest -a $BED -b galGal_ncbi_approxTSS.bed -D "b" -t "all" > $BED.annot.closest_TSS_ncbi	
	bedtools closest -a $BED -b galGal_ens_approxTSS.bed -D "b" -t "all" > $BED.annot.closest_TSS_ens
done


