## CODE TO PARSE PHYLOP OUTPUT ##
## UPDATED OCT 2018 FOR MANUSCRIPT REVISIONS ##

library(tidyverse)
setwd("~/Projects/birds/ratite_compgen/data/2018-09-FINAL/CNEE/cnee_results/phyloP_output/")

allRatite <- read.table("allRatite.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
tinamou <- read.table("Tinamou.out_neut_ver3.results", header=T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
anoDid <- read.table("anoDid.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
casuar <- read.table("Casuar.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
kiwi <- read.table("Kiwi.out_neut_ver3.results", header =T , stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
rhea <- read.table("Rhea.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
strCam <- read.table("strCam.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))

cnee_phyloP <- inner_join(allRatite, anoDid, by=c("name" = "name"), suffix=c("", ".ti")) %>%
  inner_join(anoDid, by=c("name" = "name"), suffix=c("", ".mo")) %>%
  inner_join(casuar, by=c("name" = "name"), suffix=c("", ".cd")) %>%
  inner_join(kiwi, by=c("name" = "name"), suffix=c("", ".ki")) %>%
  inner_join(rhea, by=c("name" = "name"), suffix=c("", ".rh")) %>%
  inner_join(strCam, by=c("name" = "name"), suffix=c("", ".os"))

write_tsv(cnee_phyloP, path="../cnee_phyloP.tsv")