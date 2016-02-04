#code to generate new expression estimates for chicken RNA-seq
#run with kallisto version 0.42.4 (https://pachterlab.github.io/kallisto/download.html)

#get data:
wget ftp://ftp.ncbi.nih.gov/genomes/Gallus_gallus/ARCHIVE/ANNOTATION_RELEASE.102/RNA/rna.fa.gz
gunzip rna.fa.gz
mv rna.fa galGal_NCBI.fa

#ensembl
wget ftp://ftp.ensembl.org/pub/release-83/fasta/gallus_gallus/cdna/Gallus_gallus.Galgal4.cdna.all.fa.gz
gunzip Gallus_gallus.Galgal4.cdna.all.fa.gz
mv Gallus_gallus.Galgal4.cdna.all.fa galGal_ens.fa

#make kallisto indexes
kallisto index -i galGal_ens galGal_ens.fa 
kallisto index -i galGal_NCBI galGal_NCBI.fa  

#calculate expression for each library against each index (<5 minutes per library)
for SAMPLE in HH18 HH20 HH22
do
	bzip2 -d raw_data/$SAMPLE/*.bz2
	kallisto quant -i galGal_ens -o galGal_ens_$SAMPLE --plaintext --t 4 --bias raw_data/$SAMPLE/*.fastq
	kallisto quant -i galGal_NCBI -o galGal_NCBI_$SAMPLE --plaintext --t 4 --bias raw_data/$SAMPLE/*.fastq
done
