library(dplyr)
library(tidyr)
setwd("~/Dropbox/Doctorate/RatiteManuscript/Science_Reviews/CNEE/NEW_VENN_OVERLAPS/")

chip <- read.table("~/Dropbox/Doctorate/RatiteManuscript/Science_Reviews/CNEE/NEW_VENN_OVERLAPS/update_ChIP_annotated.txt", header=T, stringsAsFactors = F) %>% tbl_df
atac <- read.table("~/Dropbox/Doctorate/RatiteManuscript/Science_Reviews/CNEE/NEW_VENN_OVERLAPS/update_ATAC_annotated.txt", header=T, stringsAsFactors = F) %>% tbl_df
rar <- read.table("~/Dropbox/Doctorate/RatiteManuscript/Science_Reviews/CNEE/NEW_VENN_OVERLAPS/update_RAR_annotated.txt", header=T, stringsAsFactors = F) %>% tbl_df


#CONVERGENT ACCELERATED CNEEs

#number of total accel cnees: 527
nrow(rar)

#number of accel cnees in only FL peaks: 21
nrow(rar %>% filter(ATAC>0,ChIP==0))

#number of accel cnees in only ChIP peaks: 67
nrow(rar %>% filter(ATAC==0,ChIP>0))

#number of accel cnees in both FL and ChIP peaks: 54
nrow(rar %>% filter(ATAC>0,ChIP>0))

#output this list of 54 for the paper
#write.table(x = rar %>% filter(ATAC>0,ChIP>0) %>% select(1:3),file = "54_RARs_Venn_Center.txt",col.names = F,row.names = F,quote = F,sep = "\t")

#number of accel cnees in neither FL and ChIP peaks: 385
nrow(rar %>% filter(ATAC==0,ChIP==0))


#FL ATAC-seq PEAKS
#number of total FL atac-seq peaks:78007
nrow(atac)

#number of FL atac-seq peaks only in cnees: 18
nrow(atac %>% filter(CNEE>0,PEAK==0))

#number of atac only in ChIP: 44383
nrow(atac %>% filter(CNEE==0,PEAK>0))

#number of atac in neither: 33550
nrow(atac %>% filter(CNEE==0,PEAK==0))

#number of atac in both: 53
nrow(atac %>% filter(CNEE>0,CNEE.PEAK>0))


#"ACTIVE" ChIP PEAKS
#number of total ChIP peaks: 102812
nrow(chip)

#number of ChIP peaks only in cnees: 36
nrow(chip %>% filter(CNEE>0,PEAK==0))

#number of ChIP only in ATAC: 40395
nrow(chip %>% filter(CNEE==0,PEAK>0))

#number of ChIP in neither: 62300
nrow(chip %>% filter(CNEE==0,PEAK==0))

#number of ChIP in both: 53
nrow(chip %>% filter(CNEE.PEAK>0,PEAK>0))
