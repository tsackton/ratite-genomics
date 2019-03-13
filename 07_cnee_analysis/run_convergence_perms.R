#this code runs the permutations to test whether convergence is greater than expected
library(tidyverse)
library(parallel)
library(ape)

#setwd("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/")

args <- commandArgs(trailingOnly = TRUE)

OUTNUM <- args[1]
CORES <- args[2]
NPERM <- args[3]
path_to_data <-args[4]

## LOAD DATA ##

cnee_orig <- read_tsv(paste0(path_to_data, "/", "final_original_cnee.tsv.gz"))
cnee_red <- read_tsv(paste0(path_to_data, "/", "final_reduced_cnee.tsv.gz"))
cnee_ext <- read_tsv(paste0(path_to_data, "/", "final_extended_cnee.tsv.gz"))

phy_orig <- read.tree(paste0(path_to_data, "/", "original.phy"))
phy_red <- read.tree(paste0(path_to_data, "/", "reduced.phy"))
phy_ext <- read.tree(paste0(path_to_data, "/", "extended.phy"))

## FUNCTIONS ##

#functions
tto <- function(x) {
  ifelse(x > 1, 1, x)
}

ctb <- function(x, cutoff=0.90, lower=FALSE) {
  if (lower) {
    ifelse(x <= cutoff, 1, 0)
  }
  else {
    ifelse(x >= cutoff, 1, 0)
  }
}

get_max <- function(x, cutoff = 0.90) {
  ifelse(max(x) > cutoff, 1, 0)
}

count_sister_taxa <- function(tree, tips) {
  #set up pairs
  pairs<-combn(tips, 2, simplify = FALSE)
  lapply(pairs, is.monophyletic, phy=tree) %>% unlist %>% sum
}

conv_sample <- function(perm=1, DF, number, tips, phy, cutoff = 0.90) {
  num_sister = 10
  while (num_sister > 0) {
    targets<-sort(sample(tips, number))
    num_sister = count_sister_taxa(phy, targets)  
  }
  #targets now has N (= number) random non-sister tips
  
  count_acc_targets <- function(tip, prob_acc, cutoff, targets) {
    targets_selected <- prob_acc[tip %in% targets]
    return(sum(targets_selected > cutoff))
  }
  
  #compute 
  
  count<-DF %>% dplyr::select(cnee, one_of(tips)) %>% 
    gather(key = "tip", value = "prob_acc", -cnee) %>%
    group_by(cnee) %>%
    summarize(total_acc = sum(prob_acc > cutoff), target_acc = count_acc_targets(tip, prob_acc, cutoff, targets)) %>%
    dplyr::filter(target_acc == number & target_acc == total_acc) %>%
    tally() %>% pull(n)
  
  return(tibble(count=count, tips=paste0(targets, collapse="-")))
}

conv_sample_cross <- function(perm=1, DF, number, tips, cross_tips, phy, cutoff = 0.90) {
  #bad progamming practice copying functions from above...oh well
  if (number > 1) {
    num_sister = 10
    while (num_sister > 0) {
      targets<-sort(sample(tips, number))
      num_sister = count_sister_taxa(phy, targets)  
    } 
  } else {
      targets<-sort(sample(tips, 1))
  }
  #targets now has N (= number) random non-sister tips
  
  count_acc_targets <- function(tip, prob_acc, cutoff, targets) {
    targets_selected <- prob_acc[tip %in% targets]
    return(sum(targets_selected > cutoff))
  }
  
  #cross target is the single "extra" species to compare
  cross_target<-sample(cross_tips, 1)
  
  #compute 
  
  count<-DF %>% dplyr::select(cnee, one_of(tips), cross_target) %>% 
    gather(key = "tip", value = "prob_acc", -cnee) %>%
    group_by(cnee) %>%
    summarize(total_acc = sum(prob_acc > cutoff), 
              target_acc = count_acc_targets(tip, prob_acc, cutoff, targets), 
              cross_target_acc = count_acc_targets(tip, prob_acc, cutoff, cross_target)) %>%
    dplyr::filter(target_acc == number, target_acc == (total_acc-1), cross_target_acc == 1) %>%
    tally() %>% pull(n)
  
  return(tibble(count=count, tips=paste0(c(targets, cross_target), collapse="-")))  
}


#EACH DATASET NEEDS DIFFERENT PROCESSING SO NEED TO DO IN SEQUENCE

##ORIGINAL (+moa, -corm)

orig_perms<-list()
red_perms<-list()
ext_perms<-list()
neo_tips_orig <- phy_orig$tip.label[1:23]
neo_tips_ext <- phy_ext$tip.label[c(1:14,16:27)]
tin_tips <- phy_ext$tip.label[35:38]

orig_perms[[1]] <- cnee_orig %>% filter(version=="gain_gap") %>% 
  mclapply(1:NPERM, conv_sample, DF=., number=3, tips=neo_tips_orig, phy=phy_orig, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain_gap", dataset="original", test="neo_conv_3") %>% 
  distinct(tips, .keep_all=TRUE)

orig_perms[[2]] <- cnee_orig %>% filter(version=="gain") %>% 
  mclapply(1:NPERM, conv_sample, DF=., number=3, tips=neo_tips_orig, phy=phy_orig, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain", dataset="original", test="neo_conv_3") %>% 
  distinct(tips, .keep_all=TRUE)

orig_perms[[3]] <- cnee_orig %>% filter(version=="gain_gap") %>% 
  mclapply(1:NPERM, conv_sample, DF=., number=2, tips=neo_tips_orig, phy=phy_orig, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain_gap", dataset="original", test="neo_conv_2") %>% 
  distinct(tips, .keep_all=TRUE)

orig_perms[[4]] <- cnee_orig %>% filter(version=="gain") %>% 
  mclapply(1:NPERM, conv_sample, DF=., number=2, tips=neo_tips_orig, phy=phy_orig, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain", dataset="original", test="neo_conv_2") %>% 
  distinct(tips, .keep_all=TRUE)

red_perms[[1]] <- cnee_red %>% filter(version=="gain_gap") %>% 
  mclapply(1:NPERM, conv_sample, DF=., number=3, tips=neo_tips_orig, phy=phy_red, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain_gap", dataset="original", test="neo_conv_3") %>% 
  distinct(tips, .keep_all=TRUE)

red_perms[[2]] <- cnee_red %>% filter(version=="gain") %>% 
  mclapply(1:NPERM, conv_sample, DF=., number=3, tips=neo_tips_orig, phy=phy_red, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain", dataset="original", test="neo_conv_3") %>% 
  distinct(tips, .keep_all=TRUE)

red_perms[[3]] <- cnee_red %>% filter(version=="gain_gap") %>% 
  mclapply(1:NPERM, conv_sample, DF=., number=2, tips=neo_tips_orig, phy=phy_red, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain_gap", dataset="original", test="neo_conv_2") %>% 
  distinct(tips, .keep_all=TRUE)

red_perms[[4]] <- cnee_red %>% filter(version=="gain") %>% 
  mclapply(1:NPERM, conv_sample, DF=., number=2, tips=neo_tips_orig, phy=phy_red, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain", dataset="original", test="neo_conv_2") %>% 
  distinct(tips, .keep_all=TRUE)

ext_perms[[1]] <- cnee_ext %>% filter(version=="gain_gap") %>% 
  mclapply(1:NPERM, conv_sample, DF=., number=3, tips=neo_tips_ext, phy=phy_ext, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain_gap", dataset="extended", test="neo_conv_3") %>% 
  distinct(tips, .keep_all=TRUE)

ext_perms[[2]] <- cnee_ext %>% filter(version=="gain") %>% 
  mclapply(1:NPERM, conv_sample, DF=., number=3, tips=neo_tips_ext, phy=phy_ext, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain", dataset="extended", test="neo_conv_3") %>% 
  distinct(tips, .keep_all=TRUE)

ext_perms[[3]] <- cnee_ext %>% filter(version=="gain_gap") %>% 
  mclapply(1:NPERM, conv_sample, DF=., number=2, tips=neo_tips_ext, phy=phy_ext, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain_gap", dataset="extended", test="neo_conv_2") %>% 
  distinct(tips, .keep_all=TRUE)

ext_perms[[4]] <- cnee_ext %>% filter(version=="gain") %>% 
  mclapply(1:NPERM, conv_sample, DF=., number=2, tips=neo_tips_ext, phy=phy_ext, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain", dataset="extended", test="neo_conv_2") %>% 
  distinct(tips, .keep_all=TRUE)

ext_perms[[5]] <- cnee_ext %>% filter(version=="gain") %>% 
  mclapply(1:NPERM, conv_sample_cross, DF=., number=1, tips=neo_tips_ext, phy=phy_ext, cross_tips = tin_tips, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain", dataset="extended", test="neo_cross_1") %>% 
  distinct(tips, .keep_all=TRUE)

ext_perms[[6]] <- cnee_ext %>% filter(version=="gain") %>% 
  mclapply(1:NPERM, conv_sample_cross, DF=., number=2, tips=neo_tips_ext, phy=phy_ext, cross_tips = tin_tips, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain", dataset="extended", test="neo_cross_2") %>% 
  distinct(tips, .keep_all=TRUE)

ext_perms[[7]] <- cnee_ext %>% filter(version=="gain") %>% 
  mclapply(1:NPERM, conv_sample_cross, DF=., number=3, tips=neo_tips_ext, phy=phy_ext, cross_tips = tin_tips, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain", dataset="extended", test="neo_cross_3") %>% 
  distinct(tips, .keep_all=TRUE)

ext_perms[[8]] <- cnee_ext %>% filter(version=="gain_gap") %>% 
  mclapply(1:NPERM, conv_sample_cross, DF=., number=1, tips=neo_tips_ext, phy=phy_ext, cross_tips = tin_tips, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain_gap", dataset="extended", test="neo_cross_1") %>% 
  distinct(tips, .keep_all=TRUE)

ext_perms[[9]] <- cnee_ext %>% filter(version=="gain_gap") %>% 
  mclapply(1:NPERM, conv_sample_cross, DF=., number=2, tips=neo_tips_ext, phy=phy_ext, cross_tips = tin_tips, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain_gap", dataset="extended", test="neo_cross_2") %>% 
  distinct(tips, .keep_all=TRUE)

ext_perms[[10]] <- cnee_ext %>% filter(version=="gain_gap") %>% 
  mclapply(1:NPERM, conv_sample_cross, DF=., number=3, tips=neo_tips_ext, phy=phy_ext, cross_tips = tin_tips, mc.preschedule = TRUE, mc.cores = CORES) %>%
  bind_rows(.id="perm") %>%
  mutate(version = "gain_gap", dataset="extended", test="neo_cross_3") %>% 
  distinct(tips, .keep_all=TRUE)

orig_perms %>% bind_rows() %>% write_tsv(paste0(path_to_data, "/", "convperms/original_perms_run", args[1], ".tsv"))
red_perms %>% bind_rows() %>% write_tsv(paste0(path_to_data, "/", "convperms/reduced_perms_run", args[1], ".tsv"))
ext_perms %>% bind_rows() %>% write_tsv(paste0(path_to_data, "/", "convperms/extended_perms_run", args[1], ".tsv"))