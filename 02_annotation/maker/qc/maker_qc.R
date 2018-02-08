setwd("~/Projects/birds/ratite_compgen/data/qc_data/")

library(tidyverse)
library(data.table)

venog<-fread("gzip -cd venog_scan.out.gz | grep -v '^#'")
venog_genes <- venog %>% mutate(gene = sub("-\\w+", "",V3)) %>% distinct(gene) %>% mutate(species = substr(gene, 1,4))

venog_genes %>% group_by(species) %>% summarize(spcount = n())

#blastp
setwd("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/blastp")
blastp<-read.table("all_besthit.out", header=F, stringsAsFactors = F)
blastp <- tbl_df(blastp)
blastp <- blastp %>% 
  rename(file = V1, query_id = V2, subj_id = V3, perc_ident = V4, length = V5, mismatch = V6, gapopen = V7, query_start = V8, query_end = V9, subj_start = V10, subj_end = V11, eval = V12, bitscore = V13) %>%
  mutate(species = sub("_galGal_blastp.out", "", file, fixed=T))

blastp %>% group_by(species) %>% summarize(count=round(n()/15732*100,2)) %>% print.data.frame()
