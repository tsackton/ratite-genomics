## CODE TO PARSE PHYLOP OUTPUT ##
## UPDATED OCT 2018 FOR MANUSCRIPT REVISIONS ##

library(tidyverse)
path_to_data<-"/Users/tim/Projects/birds/ratite_compgen/DRYAD/07_cnees/results/orig_v1_phyloP"

allRatite <- read.table(paste0(path_to_data, "/", "allRatite.out_neut_ver3.results"), header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
tinamou <- read.table(paste0(path_to_data, "/", "Tinamou.out_neut_ver3.results"), header=T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
anoDid <- read.table(paste0(path_to_data, "/", "anoDid.out_neut_ver3.results"), header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
casuar <- read.table(paste0(path_to_data, "/", "Casuar.out_neut_ver3.results"), header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
kiwi <- read.table(paste0(path_to_data, "/", "Kiwi.out_neut_ver3.results"), header =T , stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
rhea <- read.table(paste0(path_to_data, "/", "Rhea.out_neut_ver3.results"), header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
strCam <- read.table(paste0(path_to_data, "/", "strCam.out_neut_ver3.results"), header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))

cnee_phyloP <- inner_join(allRatite, anoDid, by=c("name" = "name"), suffix=c("", ".ti")) %>%
  inner_join(anoDid, by=c("name" = "name"), suffix=c("", ".mo")) %>%
  inner_join(casuar, by=c("name" = "name"), suffix=c("", ".cd")) %>%
  inner_join(kiwi, by=c("name" = "name"), suffix=c("", ".ki")) %>%
  inner_join(rhea, by=c("name" = "name"), suffix=c("", ".rh")) %>%
  inner_join(strCam, by=c("name" = "name"), suffix=c("", ".os"))

write_tsv(cnee_phyloP, path="cnee_phyloP.tsv")