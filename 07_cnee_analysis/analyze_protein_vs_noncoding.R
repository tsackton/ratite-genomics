#compare BS-REL and phyloAcc results
library(tidyverse)
library(ape)
library(jaccard)

ctb <- function(x, cutoff=0.90, lower=FALSE) {
  if (lower) {
    ifelse(x <= cutoff, 1, 0)
  }
  else {
    ifelse(x >= cutoff, 1, 0)
  }
}


#load protein matrix
prot <- read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/06_protein_coding_analysis/hyphy_bsrel/bsrel_matrix.tsv")

cnee_raw <- read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/final_extended_cnee.tsv")
cnee <- cnee_raw %>% filter(version == "gain_gap") %>% select(taeGut:strCam) %>% mutate_all(ctb)

#get tips
prot_tips <- names(prot)[1:44]
cnee_tips <- names(cnee)

#joint tips
joint_tips <- intersect(prot_tips, cnee_tips)

ratite_tips <- c("strCam", "droNov", "casCas", "anoDid", "rheAme", "rhePen", "aptHaa", "aptRow", "aptOwe")
tin_tips <- c("eudEle", "tinGut", "cryCin", "notPer")
neo_tips <- joint_tips[!(joint_tips %in% c(ratite_tips, tin_tips))]

#make tree
joint_tree <- read.tree("~/Projects/birds/ratite_compgen/ratite-genomics/06_protein_coding_analysis/final_tree_ext_proteins.nwk") %>% keep.tip(joint_tips)

#compute jaccard

compute_jaccard <- function(tip1, tip2, df, center=TRUE) {
  jaccard(df[[tip1]], df[[tip2]], center=center)
}

nt <- length(joint_tips)

cnee_sim <- tibble(tip1 = rep(joint_tips, times=nt), tip2 = rep(joint_tips, each=nt)) %>% rowwise() %>% mutate(jaccard = compute_jaccard(tip1, tip2, cnee, center=FALSE))
prot_sim <- tibble(tip1 = rep(joint_tips, times=nt), tip2 = rep(joint_tips, each=nt)) %>% rowwise() %>% mutate(jaccard = compute_jaccard(tip1, tip2, prot, center=FALSE))

all_sim <- inner_join(cnee_sim, prot_sim, by=c("tip1" = "tip1", "tip2" = "tip2"), suffix=c(".cnee", ".prot")) %>%
  distinct(jaccard.cnee, jaccard.prot, .keep_all = TRUE) %>% filter(tip1 != tip2) %>%
  mutate(pair_type = case_when(
    tip1 %in% ratite_tips & tip2 %in% ratite_tips ~ "ratite",
    tip1 %in% ratite_tips & tip2 %in% tin_tips ~ "palaeo",
    tip2 %in% ratite_tips & tip1 %in% tin_tips ~ "palaeo",
    tip1 %in% ratite_tips & tip2 %in% neo_tips ~ "neo",
    tip2 %in% ratite_tips & tip1 %in% neo_tips ~ "neo",
    TRUE ~ "other"
  ))

#trying heatmap
cnee_sim %>% ggplot(aes(tip1, tip2)) + geom_tile(aes(fill=jaccard),color="white") + scale_fill_gradient(low="white", high="steelblue") + theme(axis.text.x  = element_text(angle = 90))

prot_sim %>% ggplot(aes(tip1, tip2)) + geom_tile(aes(fill=jaccard),color="white") + scale_fill_gradient(low="white", high="steelblue") + theme(axis.text.x  = element_text(angle = 90))

#boxplot with sister taxa removed
all_sim %>% mutate(sister = is.monophyletic(joint_tree,c(tip1,tip2))) %>% 
  filter(sister == FALSE, pair_type != "other") %>%
  ggplot(aes(pair_type, jaccard.cnee)) + geom_jitter(width=0.1) 


