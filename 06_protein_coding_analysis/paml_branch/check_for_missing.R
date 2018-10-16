library(tidyverse)

#check for missing runs

setwd("~/Projects/birds/ratite_compgen/ratite-genomics/06_protein_coding_analysis/paml_branch/")

old_run <- read_tsv("aa_trees-2017-11-20.out", col_names=c("tree", "hog", "aa"), col_types = "ccc") %>% 
  filter(!is.na(hog), !is.na(aa))
new_run <- read_tsv("aa_trees-2018-10-15.out", col_names=c("tree", "hog", "aa"), col_types = "ccc") %>% 
  filter(!is.na(hog), !is.na(aa))

#get full list I should have from ancrec

full_list <- read_tsv("../paml_ancrec/paml_M0_parsed.txt.gz", col_types = "cccc?????") %>% 
  filter(!is.na(hog), model == "ancrec") %>%
  mutate(tree = paste0("tree", treenum)) %>%
  select(hog, tree, species_tree)

#merge to get missing from old run

old_run_merge <- full_list %>% left_join(old_run, by=c("tree" = "tree", "hog" = "hog")) %>% select(hog, tree, species_tree, aa)

#missing

old_run_merge %>% filter(is.na(aa)) %>% select(hog) %>% write_tsv(path="old_run_missing", col_names = FALSE)

#merge to get missing from new run

new_run_merge <- full_list %>% left_join(new_run, by=c("tree" = "tree", "hog" = "hog")) %>% select(hog, tree, species_tree, aa)

#missing

new_run_merge %>% group_by(hog) %>% summarize(good_count = sum(!is.na(aa))) %>% filter(good_count == 0) %>% select(hog) %>%
  write_tsv(path="new_run_missing", col_names = FALSE)

## RELAX ##

#Ratite RELAX_Ratite 0078 9078 tree1 all output 1.015100995437099 

#relax <- read_delim("../hyphy_relax/relax_K.all", delim = " ", col_names = c("group", "run", "dir", "hog", "tree", "model", "skip", "K", "skip2")) %>% mutate(hog = as.character(hog)) %>%
#  select(hog, tree, run, model)

#relax_merge <- full_list %>% left_join(relax, by=c("tree" = "tree", "hog" = "hog")) %>% filter(run != "RELAX_RND")

#just look at species tree

#relax_merge %>% filter(species_tree == "True") %>% group_by(hog) %>% count() %>% filter(n != 4) %>% print.data.frame

#full_list %>% count(species_tree)

#RELAX looks fine


#bs-rel

bsrel_orig <- read_tsv("../hyphy_bsrel/bsrel_res_parsed_ratites_2017-11-20.clean.txt", 
                       col_names = c("class", "tree", "hog", "tsel.s", "nsel.s", "tsel.n", "nsel.n"), col_types="ccciiii") %>% 
  filter(!is.na(hog))


bsrel_new <- read_tsv("../hyphy_bsrel/bsrel_res_parsed_ratites_2018-10-15.clean.txt", 
                       col_names = c("class", "tree", "hog", "tsel.s", "nsel.s", "tsel.n", "nsel.n"), col_types="ccciiii") %>% 
  filter(!is.na(hog))

#get missing

bso_merge <- full_list %>% left_join(bsrel_orig, by=c("tree" = "tree", "hog" = "hog"))
bsn_merge <- full_list %>% left_join(bsrel_new, by=c("tree" = "tree", "hog" = "hog"))

bso_merge %>% filter(is.na(tsel.s)) %>% select(hog) %>% write_tsv(path="../hyphy_bsrel/missing_hogs_run1", col_names = FALSE)
bsn_merge %>% filter(is.na(tsel.s)) %>% filter(tree == "tree1") %>% select(hog) %>% 
  write_tsv(path="../hyphy_bsrel/missing_hogs_run2_tree1", col_names = FALSE)
bsn_merge %>% filter(is.na(tsel.s)) %>% filter(tree == "tree2") %>% select(hog) %>% 
  write_tsv(path="../hyphy_bsrel/missing_hogs_run2_tree2", col_names = FALSE)

#how many hogs have neither tree1 nor tree2
bsn_merge %>% group_by(hog) %>% summarize(good_count = sum(!is.na(tsel.s))) %>% filter(good_count == 0) %>% select(hog) %>%
  write_tsv(path="../hyphy_bsrel/missing_hogs_run2", col_names=FALSE)

