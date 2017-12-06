#this code runs the permutations to test for GO, gene, and spatial enrichment in RARs and cRARs
library(tidyverse)
library(clusterProfiler)
library(org.Gg.eg.db)
library(parallel)
library(rlist)

#set working directory
#setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/")

#set cores and perms for whole run, can be individually modified later in code
CORES <- 32
PERMS <- 12000

#load data
cnee <- read_tsv("cnees.tsv")

### GO ENRICHMENT HERE ###
#GO enrichment - four tests: accel .1s, accel .1s & conv .1s, accel .2s, accel. 1s & conv .1s & ratite_loss_cons_min.mat >= 2
#this code gets the real results
background <- cnee %>% filter(gene != ".") %>% dplyr::select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi)
set1 <- cnee %>% filter(gene != ".", ratite_accel.1, ratite_spec.1) %>% dplyr::select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi)
set2 <- cnee %>% filter(gene != ".", ratite_accel.2, ratite_spec.2) %>% dplyr::select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi)
set3 <- cnee %>% filter(gene != ".", ratite_accel.1, ratite_spec.1, ratite_conv.1) %>% dplyr::select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi)
set4 <- cnee %>% filter(gene != ".", ratite_accel.1, ratite_spec.1, ratite_conv.1, ratite_loss_cons_min.mat >= 2) %>% dplyr::select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi)
inputs <- list("set1" = set1, "set2" = set2, "set3" = set3, "set4" = set4)
calc_enrich <- function(targetset, background,ont) { enrichGO(targetset$ncbi,'org.Gg.eg.db',pvalueCutoff=1.5,qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, pAdjustMethod="none",universe=background,keytype="ENTREZID",ont=ont) }
bp_res_real <- lapply(inputs, calc_enrich, background=background$ncbi, ont="BP") %>% 
  lapply(slot, name="result") %>% 
  dplyr::bind_rows(.id = "set")
mf_res_real <-  lapply(inputs, calc_enrich, background=background$ncbi, ont="MF") %>% 
  lapply(slot, name="result") %>% 
  dplyr::bind_rows(.id = "set")
merged_mf_terms <- mf_res_real %>% dplyr::distinct(ID)
merged_bp_terms <- bp_res_real %>% dplyr::distinct(ID)

#get counts of CNEEs in each set
input_counts<-list("set1" = cnee %>% filter(gene != ".", ratite_accel.1, ratite_spec.1) %>% count %>% pull(n),
                   "set2" = cnee %>% filter(gene != ".", ratite_accel.2, ratite_spec.2) %>% count %>% pull(n),
                   "set3" = cnee %>% filter(gene != ".", ratite_accel.1, ratite_spec.1, ratite_conv.1) %>% count %>% pull(n),
                   "set4" = cnee %>% filter(gene != ".", ratite_accel.1, ratite_spec.1, ratite_conv.1, ratite_loss_cons_min.mat >= 2) %>% count  %>% pull(n))

#for each GO term in the merged go list, want to compute permutations: input is the merged term list and the counts

get_go_perm <- function(DF, samples, golist, ont) {
  rand <- DF %>% sample_n(samples) %>% filter(gene != ".") %>% 
    dplyr::select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
    distinct(ncbi)
  background <- DF %>% filter(gene != ".") %>% dplyr::select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi)
  rand_go_bp <- calc_enrich(targetset=rand, background=background$ncbi, ont=ont)
  golist %>% left_join(rand_go_bp@result, by=c("ID" = "ID")) %>% separate(GeneRatio, into=c("target_in", "target_total")) %>% 
    separate(BgRatio, into=c("bg_in", "bg_total")) %>%
    mutate(newpval = ifelse(is.na(pvalue), 1, pvalue), logp.perm = -log10(newpval)) %>% 
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total), bg_frac = as.numeric(bg_in)/as.numeric(bg_total)) %>%
    dplyr::select(ID, logp.perm, target_frac, bg_frac) %>% arrange(ID)
}

#to do one permutation of the full set
get_one_perm_set <- function(perm, input, DF, golist, ont) {
  try(lapply(input, get_go_perm, DF=DF, golist=golist, ont=ont) %>%
  dplyr::bind_rows(.id="set"), TRUE)
}

#do permutations in parallel
para_cores <- CORES
num_perms <- PERMS

perm_bp <- mclapply(1:num_perms, get_one_perm_set, input=input_counts, DF=cnee, golist=merged_bp_terms, ont="BP", mc.cores=para_cores, mc.preschedule = FALSE) %>%
  list.filter(class(.) == "data.frame") %>% 
  dplyr::bind_rows(.id="perm")
perm_mf <- mclapply(1:num_perms, get_one_perm_set, input=input_counts, DF=cnee, golist=merged_mf_terms, ont="MF", mc.cores=para_cores, mc.preschedule = FALSE) %>%
  list.filter(class(.) == "data.frame") %>% 
  dplyr::bind_rows(.id="perm")

#write out
write_tsv(perm_bp, path="perm_bp_results.tsv")
write_tsv(perm_mf, path="perm_mf_results.tsv")
write_tsv(bp_res_real, path="obs_bp_results.tsv")
write_tsv(mf_res_real, path="obs_mf_results.tsv")

#end GO results#

### GENE ENRICHMENT HERE ###
#for gene enrichment tests, the idea is to randomly permute each set and get count of CNEEs per gene
#strategy here is to make an indicator variable, compute "real" T/F per gene, and then shuffle indicator variable

get_gene_counts <- function(DF, indicator) {
#  indcol <- enquo(indicator)
#  indcol <- indicator
  DF %>% filter(gene != ".") %>% mutate(in_target = !!indicator) %>% count(gene, in_target) %>% filter(!is.na(in_target)) %>% spread(in_target, n, fill=0, drop=FALSE, sep="_")
}

perm_gene_counts <- function(perm, DF, indicator) {
#  indcol <- enquo(indicator)
  DF %>% filter(gene != ".") %>% mutate(rand = sample(!!indicator)) %>% count(gene, rand) %>% filter(!is.na(rand)) %>% spread(rand, n, fill=0, drop=FALSE, sep="_")
}

#make real dataset
#first add sets to cnee data frame
cnee <- cnee %>% mutate(set1 = ratite_accel.1 & ratite_spec.1) %>% 
  mutate(set2 = ratite_accel.2 & ratite_spec.2) %>% 
  mutate(set3 = ratite_accel.1 & ratite_spec.1 & ratite_conv.1) %>% 
  mutate(set4 = ratite_accel.1 & ratite_spec.1 & ratite_conv.1 & ratite_loss_cons_min.mat >= 2)

gene_counts_obs <- lapply(c(quo(set1),quo(set2),quo(set3),quo(set4)), get_gene_counts, DF=cnee) %>% bind_rows(.id="set")

para_cores <- CORES
num_perms <- PERMS

gene_counts_perm_list <- list(set1 = mclapply(1:num_perms, perm_gene_counts, DF=cnee, indicator=quo(set1), mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"),
                                   set2 = mclapply(1:num_perms, perm_gene_counts, DF=cnee, indicator=quo(set2), mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"),
                                   set3 = mclapply(1:num_perms, perm_gene_counts, DF=cnee, indicator=quo(set3), mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"),
                                   set4 = mclapply(1:num_perms, perm_gene_counts, DF=cnee, indicator=quo(set4), mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"))

gene_counts_perm <- gene_counts_perm_list %>% bind_rows(.id="set")

write_tsv(gene_counts_perm, path="perm_gene_count_results.tsv")
write_tsv(gene_counts_obs, path="obs_gene_count_results.tsv")

#end GENE PERMUTATIONS#
