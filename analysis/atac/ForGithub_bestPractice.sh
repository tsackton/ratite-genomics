#Initial processing, mapping and peak calling follows the ATAC-seq best practices from Harvard RC computing as written by John Gaspar.  
#Programs and pipeline are available here: https://github.com/harvardinformatics/ATAC-seq
#Chicken FL and HL were sequenced on a single run, but Keel, Sternum, and Pec were across two lanes.  The first step was to concatenate the raw reads from these libraries

#raw reads were found in KSP_Run1 and KSP_Run2.  The concatenated reads were sent to Oct31_BestPractices 
for FILE in $(ls KSP_Run1/C*.gz); do cat KSP_Run1/${FILE##*/} KSP_Run2/${FILE##*/} > Oct31_BestPractices/cat_${FILE##*/}; done &

#adapters were trimmed using NGmerge
for FILE in $(ls C*.R1*gz); do FILE2=$(echo $FILE | awk -F'[.]' '{print $1}'); sbatch NGmerge_20.sh $FILE2; sleep 1; done

#the files were then mapped against galGal4 using bowtie2
for temp in $(ls 20_Cat*_1*.gz); do FILE2=$(echo $temp | awk -F '[_]' '{print "20_Cat_"$3"_"$4"_2.fastq.gz"}'); sbatch bowtie2_X2000.sh $temp $FILE2 ../galGal; sleep 1; done
for temp in $(ls 20_Chick*_1*.gz); do FILE2=$(echo $temp | awk -F '[_]' '{print "20_"$2"_"$3"_2.fastq.gz"}'); sbatch bowtie2_X2000.sh $temp $FILE2 ../galGal; sleep 1; done

#upon completion, mapping quality was checked by examining the last 15 lines of each error file  
tail -n 15 *err

#picard was used to remove duplicate reads
for FILE in $(ls *bam); do sbatch picard_remove_dups.sh $FILE; sleep 1; done

#and John Gaspar's custom removeChrom was used to remove the chicken mitochondrial chromosome
for FILE in $(ls no_dups_20_C*bam); do sbatch mito_gone.sh $FILE; sleep 1; done

#on a test/interactive node, a new bam index for each file was generated at this point
srun -p test --pty --mem 2000 -t 0-02:00 /bin/bash
source new-modules.sh
module load samtools/0.1.19-fasrc01
for FILE in $(ls noMito_no_dups_20_C*bam); do samtools index $FILE; done

#peaks were called using MACS2
#first, a relaxed significance threshold (p 0.05) was utilized for each individual library (3 bio reps x 8 tissues = 24)
for FILE in $(ls noMito_no_dups_20_C*bam); do sbatch c_indy_macs_bdg_p05_dupsRetained.sh $FILE; sleep 1; done

#then, to increase the confidence in peak boundaries, a strict cutoff (q 0.05) was used and peaks were called again using all three libraries for each bio rep (1 pool x 8 = 8)
#important to note here that for keel, sternum and pec, e9 libraries are numbered 3, 7, and 8, while e10 libraries are numbered 4, 5, and 6
sbatch c_pool_macs_bdg_q05_dupsRetained.sh noMito_no_dups_20_Cat_C-4-PecMus_S12_1.fastq.gz.bam noMito_no_dups_20_Cat_C-5-PecMus_S15_1.fastq.gz.bam noMito_no_dups_20_Cat_C-6-PecMus_S18_1.fastq.gz.bam
sbatch c_pool_macs_bdg_q05_dupsRetained.sh noMito_no_dups_20_Cat_C-3-PecMus_S3_1.fastq.gz.bam noMito_no_dups_20_Cat_C-7-PecMus_S6_1.fastq.gz.bam noMito_no_dups_20_Cat_C-8-PecMus_S9_1.fastq.gz.bam
sbatch c_pool_macs_bdg_q05_dupsRetained.sh noMito_no_dups_20_Cat_C-4-Sternum_S11_1.fastq.gz.bam noMito_no_dups_20_Cat_C-5-Sternum_S14_1.fastq.gz.bam noMito_no_dups_20_Cat_C-6-Sternum_S17_1.fastq.gz.bam
sbatch c_pool_macs_bdg_q05_dupsRetained.sh noMito_no_dups_20_Cat_C-3-Sternum_S2_1.fastq.gz.bam noMito_no_dups_20_Cat_C-7-Sternum_S5_1.fastq.gz.bam noMito_no_dups_20_Cat_C-8-Sternum_S8_1.fastq.gz.bam
sbatch c_pool_macs_bdg_q05_dupsRetained.sh noMito_no_dups_20_Cat_C-4-Keel_S10_1.fastq.gz.bam noMito_no_dups_20_Cat_C-5-Keel_S13_1.fastq.gz.bam noMito_no_dups_20_Cat_C-6-Keel_S16_1.fastq.gz.bam
sbatch c_pool_macs_bdg_q05_dupsRetained.sh noMito_no_dups_20_Cat_C-3-Keel_S1_1.fastq.gz.bam noMito_no_dups_20_Cat_C-7-Keel_S4_1.fastq.gz.bam noMito_no_dups_20_Cat_C-8-Keel_S7_1.fastq.gz.bam
sbatch c_pool_macs_bdg_q05_dupsRetained.sh noMito_no_dups_20_C*HL*bam
sbatch c_pool_macs_bdg_q05_dupsRetained.sh noMito_no_dups_20_C*FL*bam

#the peaks from the relaxed run were then filtered to the same strict cutoff (q 0.05) used in the pooled runs:
for FILE in $(ls noM*p05*narrowPeak); do awk '($9  >= 1.30102999566)' $FILE > filtered_${FILE}; done

#on a test/interactive node, bedtools annotate was used to link our pooled peaks to the filtered peaks from each individual library
#as above we are grouping 3, 7, and 8 for e9 and 4, 5, 6, for e10
srun -p test --pty --mem 20000 -t 0-02:00 /bin/bash
source new-modules.sh
module load bedtools2/2.26.0-fasrc01
bedtools annotate -i noMito_no_dups_20_Chick1-TR2-FL_S1_1.fastq.gz.bam_macs2_pool_keepDup_q05_peaks.narrowPeak -files filtered_noMito_no_dups_20_C*FL* -names rep1,rep2,rep3 > FL_annotate.txt &
bedtools annotate -i noMito_no_dups_20_Chick1-TR2-HL_S2_1.fastq.gz.bam_macs2_pool_keepDup_q05_peaks.narrowPeak -files filtered_noMito_no_dups_20_C*HL* -names rep1,rep2,rep3 > HL_annotate.txt &
bedtools annotate -i noMito_no_dups_20_Cat_C-3-Keel_S1_1.fastq.gz.bam_macs2_pool_keepDup_q05_peaks.narrowPeak -files filtered_noMito_no_dups_20_Cat_C-3-Keel_S1_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak filtered_noMito_no_dups_20_Cat_C-7-Keel_S4_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak filtered_noMito_no_dups_20_Cat_C-8-Keel_S7_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak -names rep1,rep2,rep3 > Keel9_annotate.txt &
bedtools annotate -i noMito_no_dups_20_Cat_C-4-Keel_S10_1.fastq.gz.bam_macs2_pool_keepDup_q05_peaks.narrowPeak -files filtered_noMito_no_dups_20_Cat_C-4-Keel_S10_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak filtered_noMito_no_dups_20_Cat_C-5-Keel_S13_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak filtered_noMito_no_dups_20_Cat_C-6-Keel_S16_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak -names rep1,rep2,rep3 > Keel10_annotate.txt &
bedtools annotate -i noMito_no_dups_20_Cat_C-3-PecMus_S3_1.fastq.gz.bam_macs2_pool_keepDup_q05_peaks.narrowPeak -files filtered_noMito_no_dups_20_Cat_C-3-PecMus_S3_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak filtered_noMito_no_dups_20_Cat_C-7-PecMus_S6_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak filtered_noMito_no_dups_20_Cat_C-8-PecMus_S9_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak -names rep1,rep2,rep3 > Pec9_annotate.txt &
bedtools annotate -i noMito_no_dups_20_Cat_C-4-PecMus_S12_1.fastq.gz.bam_macs2_pool_keepDup_q05_peaks.narrowPeak -files filtered_noMito_no_dups_20_Cat_C-4-PecMus_S12_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak filtered_noMito_no_dups_20_Cat_C-5-PecMus_S15_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak filtered_noMito_no_dups_20_Cat_C-6-PecMus_S18_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak -names rep1,rep2,rep3 > Pec10_annotate.txt &
bedtools annotate -i noMito_no_dups_20_Cat_C-3-Sternum_S2_1.fastq.gz.bam_macs2_pool_keepDup_q05_peaks.narrowPeak -files filtered_noMito_no_dups_20_Cat_C-3-Sternum_S2_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak filtered_noMito_no_dups_20_Cat_C-7-Sternum_S5_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak filtered_noMito_no_dups_20_Cat_C-8-Sternum_S8_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak -names rep1,rep2,rep3 > Strn9_annotate.txt &
bedtools annotate -i noMito_no_dups_20_Cat_C-4-Sternum_S11_1.fastq.gz.bam_macs2_pool_keepDup_q05_peaks.narrowPeak -files filtered_noMito_no_dups_20_Cat_C-4-Sternum_S11_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak filtered_noMito_no_dups_20_Cat_C-5-Sternum_S14_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak filtered_noMito_no_dups_20_Cat_C-6-Sternum_S17_1.fastq.gz.bam_macs2_keepDup_p05_peaks.narrowPeak -names rep1,rep2,rep3 > Strn10_annotate.txt &

#bedtools annotate output provides the relative overlap of the pooled peak against each individual filtered peak.
#peaks were considered strict and carried forward in the analysis so long as the overlap was greater than 0 for all 3 biological replicates
for FILE in $(ls *_annotate.txt); do awk '($11 > 0)' $FILE | awk '($12 > 0)' - | awk '($13 > 0)' - > strict_${FILE}; done

#bedtools annotate was also utilized to link generate the file used in candidate filtering steps by annotating our CNEEs with these strict ATAC peaks and the Seki et al. (2017) ChIP-seq peaks
bedtools annotate -i final_cnees_long.bed -files sort_*bed3 St*txt -names FL_strict,HL_strict,Keel10_strict,Keel9_strict,Pec10_strict,Pec9_strict,Strn10_strict,Strn9_strict,St16_whole_H3K27ac,St16_whole_H3K27me3,St16_whole_H3K4me1,St21_limb_H3K27ac,St21_limb_H3K27me3,St21_limb_H3K4me1,St21_whole_H3K27ac,St21_whole_H3K27me3,St21_whole_H3K4me1,St32_limb_H3K27ac,St32_limb_H3K27me3,St32_limb_H3K4me1,St32_whole_H3K27ac,St32_whole_H3K27me3,St32_whole_H3K4me1 > tissue_specific_Nov9.txt