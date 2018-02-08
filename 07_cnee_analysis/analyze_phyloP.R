library(dplyr)
library(tidyr)

setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/")

#read tip tests

#make list of all targets [tips and internal] to read in
targets<-read.table("phylop_target_list.txt", header=F, stringsAsFactors = F)
results<-list()
for (target in targets$V1) {
  infile=paste0(getwd(), "/phyloP_all_branches/", target, ".phyloP.out.gz")
  temp<-read.delim(infile)
  names(temp)[1]="cnee"
  temp$sp=target
  temp$qval=p.adjust(temp$pval, method="fdr")
  temp$padj=p.adjust(temp$pval, method="holm")
  results[[target]]=temp
}

phylop.fdr<-data.frame(cnee=results$anaPla$cnee)
for (target in targets$V1) {
  phylop.fdr[,target] = as.numeric(results[[target]]$qval<0.05)
}
row.names(phylop.fdr)<-phylop.fdr$cnee

phylop.holm<-data.frame(cnee=results$anaPla$cnee)
for (target in targets$V1) {
  phylop.holm[,target] = as.numeric(results[[target]]$padj<0.05)
}
row.names(phylop.holm)<-phylop.holm$cnee

#read clade tests
setwd("~/Projects/birds/ratite_compgen/data/accelTests/withMoa/")
allRatite <- read.table("allRatite.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
tinamou <- read.table("Tinamou.out_neut_ver3.results", header=T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
anoDid <- read.table("anoDid.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
casuar <- read.table("Casuar.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
kiwi <- read.table("Kiwi.out_neut_ver3.results", header =T , stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
rhea <- read.table("Rhea.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
strCam <- read.table("strCam.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
