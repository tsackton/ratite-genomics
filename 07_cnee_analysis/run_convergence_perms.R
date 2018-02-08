#this code runs the permutations to test whether convergence is greater than expected
library(tidyverse)
library(parallel)
library(ape)
#
setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/")
#
#we'll work in trios as that is conservative, using either the .mat or .prob
#for ratites, use neognaths as control

CORES <- 32
NPERM <- 20000

cnee <- read_tsv("cnees.tsv")

#first set up neognath tip names
neo_mat_tips <- cnee %>% dplyr::select(ficAlb.mat:anaPla.mat) %>% names
neo_prob_tips <- cnee %>% dplyr::select(ficAlb.prob:anaPla.prob) %>% names

#add total number of losses
cnee$neo_loss_mat_full <- cnee %>% dplyr::select(ficAlb.mat:anaPla.mat) %>% rowSums(.)
cnee$neo_loss_prob_full <- cnee %>% dplyr::select(ficAlb.prob:anaPla.prob) %>% rowSums(.)

#read tree
phy<-read.tree(file="final_neut_tree.nwk")

count_sister_taxa <- function(tree, tips) {
  #set up pairs
  pairs<-combn(tips, 2, simplify = FALSE)
  lapply(pairs, is.monophyletic, phy=phy) %>% unlist %>% sum
}

#sample
#wrap in function
conv_sample <- function(perm=1, DF, number, tips) {
  num_sister = 10
  while (num_sister > 0) {
    targets<-sort(sample(tips, number))
    tip_names<-sub(".mat", "", targets, fixed=TRUE)
    num_sister = count_sister_taxa(phy, tip_names)  
  }
  count<-DF %>% dplyr::select(one_of(targets),neo_loss.mat,neo_loss_mat_full) %>% 
    mutate(selected_loss = rowSums(.[,1:3])) %>% 
    dplyr::filter(selected_loss == 3 & neo_loss_mat_full == selected_loss) %>%
    tally() %>% pull(n)
  return(tibble(count=count, tips=paste0(targets, collapse="-")))
}

#perm_conv <- lapply(1:100, conv_sample, DF=cnee, number=3, tips=neo_mat_tips) %>% bind_rows(.id="perm")
perm_conv <- mclapply(1:NPERM, conv_sample, DF=cnee, number=3, tips=neo_mat_tips, mc.preschedule = TRUE, mc.cores = CORES) %>% bind_rows(.id="perm")
perm_conv_clean <- perm_conv %>% distinct(tips, .keep_all=TRUE)

write_tsv(perm_conv_clean, path="perm_convergence.tsv")

#get real ratite results
get_ratite_count <- function(DF, species) {
  target = quo(paste0(species, ".mat", collapse=""))
  DF %>% filter(neo_loss.mat == 0, tinGut.mat == 0, notPer.mat==0, eudEle.mat==0, cryCin.mat==0) %>% dplyr::select(strCam.mat, anoDid.mat, !!target) %>% 
    mutate(selected_loss = rowSums(.[,1:3])) %>% 
    dplyr::filter(selected_loss == 3) %>%
    tally() %>% pull(n)
}

ratite_conv_counts <- data.frame(third_species = c("rheAme", "rhePen", "casCas", "droNov", "aptHaa", "aptOwe", "aptRow"), 
                                 conv_count = lapply(c("rheAme", "rhePen", "casCas", "droNov", "aptHaa", "aptOwe", "aptRow"), get_ratite_count, DF=cnee) %>% unlist)

write_tsv(ratite_conv_counts, path="obs_convergence.tsv")
