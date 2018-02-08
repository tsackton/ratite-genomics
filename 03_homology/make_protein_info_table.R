library(dplyr)
library(tidyr)

#update with info here: /n/holylfs/LABS/edwards_lab/Users/acloutier/oma_homology/HOG_final_alignment_seqids

setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/")

all_hogs<-read.table("alignment_summary_stats", sep="\t", header=T, comment.char="", na.strings="n/a")
all_hogs$hog=sub("HOG2_", "", all_hogs$Locus)
all_hogs<-all_hogs[,c(51,3,5,6,7,8)]
names(all_hogs)=c("hog", "align_len", "avg_gapped_len", "avg_ident", "gaps_bp_tot", "gaps_bp_avg")

hog_align_info<-read.table("good_PAML_HOG_protein_transcript_info", header=T, sep="\t")
hog_align_info$hog = sub("HOG2_", "", hog_align_info$HOG)
hog_align_info = hog_align_info[,c("hog", "Species")]

hog_sp_ct<-count(as.tbl(hog_align_info), hog, Species)
hog_sp_ct_final<-tidyr::spread(hog_sp_ct, Species, n)

hog_info<-merge(hog_sp_ct_final, all_hogs, by="hog")

species=names(hog_info)[2:40]
hog_info$missing_ct = apply(hog_info[,species], 1, function(x) sum(is.na(x)))
hog_info$dup_ct = apply(hog_info[,species], 1, function(x) sum(x[!is.na(x)] > 1))

write.table(hog_info, file="all_hog_info.tsv", sep="\t", quote=F, row.names = F)
