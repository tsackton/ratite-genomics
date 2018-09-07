#screen GFFs for gene models overlapping masked sequence
#working dir: /n/holylfs/LABS/edwards_lab/tsackton/RATITE_PAPER_DATA_FREEZE/FINAL_RELEASE/02_annotation/masked_versions
#setup default environment:
source ~/default_env.sh 

#copy masks
cp ../../01_assembly/genbank_submission/ver1_masks/*.bed .

#get all elements that overlap masked regions

for SP in aptHaa aptOwe aptRow casCas cryCin notPer eudEle droNov rheAme rhePen;
do
   bedtools intersect -wa -a ../final_gffs/$SP.genome.gff -b ${SP}_ncbi_mask.bed > $SP.masked_models.bed
done

#make list of potentially problematic genes (anything with exon in the masked_models bed):

cat *.masked_models.bed | awk '{if ($3 == "CDS") print}' | cut -f9 | perl -p -e 's/^.*Parent=(\w+)-(\w+);$/$1/' | sort | uniq > flagged_models.txt

#check flagged models against hogs (good paml transcript file):

cp ../../../../ratite-genomics/03_homology/good_PAML_HOG_protein_transcript_info .
for TEST in $(cat flagged_models.txt); do grep $TEST good_PAML_HOG_protein_transcript_info >> flagged_hogs.txt; done

