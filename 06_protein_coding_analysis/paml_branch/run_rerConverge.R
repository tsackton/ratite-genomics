#use RERconverge
library(tidyverse)
library(devtools)

#using my fork of RERconverge that fixes a few small bugs
#install_github("tsackton/RERconverge")

library(RERconverge)

setwd("~/Projects/birds/ratite_compgen/ratite-genomics/06_protein_coding_analysis/paml_branch/")

## INITIAL SET UP ##

#get hog to gene id
hog_to_gene <- read.table("../../03_homology/HOG_final_alignment_seqids", header=T, stringsAsFactors =F)
hog_to_gene <- hog_to_gene %>% tbl_df %>%
  mutate(hogid = as.integer(sub("HOG2_", "", HOG, fixed=T))) %>% 
  mutate(hogname = sub("HOG2_", "HOG", HOG, fixed = T)) %>%
  filter(Taxon == "galGal") %>%
  dplyr::select(-HOG, -Fasta_seqid, -Taxon, -Transcript) %>%
  group_by(hogid) %>%
  mutate(gene = paste0(Gene, sep="", collapse=";")) %>%
  dplyr::select(-Gene)


#make trait trees
#make species labelled tree
full_tree <- read.tree("../final_tree_ext_proteins.nwk")

#prune to remove non-bird outgroups and set branch lengths to 1 
#dropping egrGar because it is in a different place in the species tree used for aaml
ext_tree <- drop.tip(full_tree, c("allMis", "chrPic", "anoCar", "egrGar")) %>% compute.brlen(0)
orig_tree <- drop.tip(full_tree, c("allMis", "chrPic", "anoCar", "egrGar", "nanAur", "nanBra", "nanHar", "uriPel")) %>% compute.brlen(0)
red_tree <- drop.tip(full_tree, c("allMis", "chrPic", "anoCar", "egrGar", "nanAur", "nanBra", "nanHar", "uriPel", "anoDid")) %>% compute.brlen(0)

### MAKE CORRECT PHENOTYPIC TREES ## 

#EXTENDED TREE -- three versions, flightless, ratite, and vocal learners

#plot
plot(ext_tree, use.edge.length = F)
edgelabels()

#make flightless tree
flightless_ext <- ext_tree
flightless_ext$edge.length[sort(c(84,83,41))] = 1
flightless_ext$edge.length[sort(c(67,68,65,66,64))] = 1
flightless_ext$edge.length[sort(c(69,70,71))] = 1
flightless_ext$edge.length[sort(c(72,73,74))] = 1
flightless_ext <- unroot(flightless_ext)

#make ratite tree
ratite_ext <- flightless_ext
ratite_ext$edge.length[40] = 0

#set edges for vocal learners to 0
vl_ext <- ext_tree
vl_ext$edge.length[sort(c(20,49))] = 1
vl_ext$edge.length[sort(c(9:19))] = 1
vl_ext <- unroot(vl_ext)


#ORIGINAL TREE -- two versions, ratite and vocal learners

#plot
plot(orig_tree, use.edge.length = F)
edgelabels()

#make ratite tree
ratite_orig <- orig_tree
ratite_orig$edge.length[sort(c(76,75))] = 1
ratite_orig$edge.length[sort(c(60:56))] = 1
ratite_orig$edge.length[sort(c(64:66))] = 1
ratite_orig$edge.length[sort(c(61:63))] = 1
ratite_orig <- unroot(ratite_orig)

#set edges for vocal learners to 0
vl_orig <- orig_tree
vl_orig$edge.length[sort(c(41,20))] = 1
vl_orig$edge.length[sort(c(9:19))] = 1
vl_orig <- unroot(vl_orig)


#REDUCED TREE -- two versions, ratite and vocal learners

#plot
plot(red_tree, use.edge.length = F)
edgelabels()

#make ratite tree
ratite_red <- red_tree
ratite_red$edge.length[sort(c(74))] = 1
ratite_red$edge.length[sort(c(60:56))] = 1
ratite_red$edge.length[sort(c(64:66))] = 1
ratite_red$edge.length[sort(c(61:63))] = 1
ratite_red <- unroot(ratite_red)

#set edges for vocal learners to 0
vl_red <- red_tree
vl_red$edge.length[sort(c(41,20))] = 1
vl_red$edge.length[sort(c(9:19))] = 1
vl_red <- unroot(vl_red)

# check all trees
plot(vl_ext)
plot(vl_orig)
plot(vl_red)

plot(ratite_ext)
plot(ratite_orig)
plot(ratite_red)

plot(flightless_ext)

#load genes - two datasets
if (file.exists("aaSpTreesRER_ext.Robj")) {
  load("aaSpTreesRER_ext.Robj")
} else {
  aaSpTreesExt <- readTrees("aa_trees-2018-09-06_RERconverge.txt")
  save(aaSpTreesExt, file="aaSpTreesRER_ext.Robj")
}

if (file.exists("aaSpTreesRER_red.Robj")) {
  load("aaSpTreesRER_red.Robj")
} else {
  aaSpTreesRed <- readTrees("aa_trees-2017-11-20_RERconverge.txt")
  save(aaSpTreesRed, file="aaSpTreesRER_red.Robj")
}

#make RER objects: extended and original from aaSpTreeExt, reduced from aaSpTreeRed

if (file.exists("RERext.Robj")) {
  load("RERext.Robj")
} else {
  RERext <- RERconverge::getAllResiduals(aaSpTreesExt, useSpecies=ext_tree$tip.label, transform = "sqrt", weighted = T, scale = T)
  save(RERext, file="RERext.Robj")
}

if (file.exists("RERorig.Robj")) {
  load("RERorig.Robj")
} else {
  RERorig <- RERconverge::getAllResiduals(aaSpTreesExt, useSpecies=orig_tree$tip.label, transform = "sqrt", weighted = T, scale = T)
  save(RERorig, file="RERorig.Robj")
}

if (file.exists("RERred.Robj")) {
  load("RERred.Robj")
} else {
  RERred <- RERconverge::getAllResiduals(aaSpTreesRed, useSpecies=red_tree$tip.label, transform = "sqrt", weighted = T, scale = T)
  save(RERred, file="RERred.Robj")
}

#reset par due to RERconverge issues
dev.off()

### DONE INITIAL SETUP ###


#make paths -- 7 total
vl_ext_paths <- tree2Paths(vl_ext, aaSpTreesExt)
vl_orig_paths <- tree2Paths(vl_orig, aaSpTreesExt)
vl_red_paths <- tree2Paths(vl_red, aaSpTreesRed)
ratite_ext_paths <- tree2Paths(ratite_ext, aaSpTreesExt)
ratite_orig_paths <- tree2Paths(ratite_orig, aaSpTreesExt)
ratite_red_paths <- tree2Paths(ratite_red, aaSpTreesRed)
flightless_ext_paths <- tree2Paths(flightless_ext, aaSpTreesExt)

### END MAKE PHENOTYPE TREES ### 

## RUN RER - 7 runs ###

res_vl_ext <- correlateWithBinaryPhenotype(RERext, vl_ext_paths, min.sp=10, min.pos=3, weighted = "auto") %>% as_tibble(rownames="gene")
res_vl_orig <- correlateWithBinaryPhenotype(RERorig, vl_orig_paths, min.sp=10, min.pos=3, weighted = "auto") %>% as_tibble(rownames="gene")
res_vl_red <- correlateWithBinaryPhenotype(RERred, vl_red_paths, min.sp=10, min.pos=3, weighted = "auto") %>% as_tibble(rownames="gene")
res_ratite_ext <- correlateWithBinaryPhenotype(RERext, ratite_ext_paths, min.sp=10, min.pos=3, weighted = "auto") %>% as_tibble(rownames="gene")
res_ratite_orig <- correlateWithBinaryPhenotype(RERorig, ratite_orig_paths, min.sp=10, min.pos=3, weighted = "auto") %>% as_tibble(rownames="gene")
res_ratite_red <- correlateWithBinaryPhenotype(RERred, ratite_red_paths, min.sp=10, min.pos=3, weighted = "auto") %>% as_tibble(rownames="gene")
res_flightless_ext <- correlateWithBinaryPhenotype(RERext, flightless_ext_paths, min.sp=10, min.pos=3, weighted = "auto") %>% as_tibble(rownames="gene")

res_RER_all <- bind_rows(list(vl_ext = res_vl_ext, vl_orig = res_vl_orig, vl_red = res_vl_red, rat_ext = res_ratite_ext, rat_orig = res_ratite_orig, rat_red = res_ratite_red, fl_ext = res_flightless_ext), .id = "run") %>% separate(run, into=c("targets", "tree")) %>% mutate(direction = ifelse(Rho < 0, "down", "up"))
write_tsv(res_RER_all, path="all_RER_results.tsv")



## END RUN RER ##

## INITIAL SUMMARIES ##

res_RER_all %>% mutate(direction = case_when(p.adj < 0.05 & Rho < 0 ~ "down", p.adj < 0.05 & Rho > 0 ~ "up", TRUE ~ "ns")) %>%
  group_by(targets, tree) %>% count(direction) %>% spread(direction, n)

res_RER_all %>% mutate(direction = case_when(p.adj < 0.01 & Rho < 0 ~ "down", p.adj < 0.01 & Rho > 0 ~ "up", TRUE ~ "ns")) %>%
  group_by(targets, tree) %>% count(direction) %>% spread(direction, n)

res_RER_all %>% mutate(direction = case_when(P < 0.05 & Rho < 0 ~ "down", P < 0.05 & Rho > 0 ~ "up", TRUE ~ "ns")) %>%
  group_by(targets, tree) %>% count(direction) %>% spread(direction, n) %>% mutate(prop_sig = (down + up) / (down + up + ns))

res_RER_all %>% mutate(direction = case_when(P < 0.01 & Rho < 0 ~ "down", P < 0.01 & Rho > 0 ~ "up", TRUE ~ "ns")) %>%
  group_by(targets, tree) %>% count(direction) %>% spread(direction, n) %>% mutate(prop_sig = (down + up) / (down + up + ns))

res_RER_all %>% filter(tree == "orig", direction=="up") %>% ggplot(aes(P, col=targets)) + geom_freqpoly(binwidth=0.05)
res_RER_all %>% filter(tree == "orig", direction=="down") %>% ggplot(aes(P, col=targets)) + geom_freqpoly(binwidth=0.05)
res_RER_all %>% filter(tree == "ext", direction=="up") %>% ggplot(aes(P, col=targets)) + geom_freqpoly(binwidth=0.05)
res_RER_all %>% filter(tree == "ext", direction=="down") %>% ggplot(aes(P, col=targets)) + geom_freqpoly(binwidth=0.05)
res_RER_all %>% filter(tree == "red", direction=="up") %>% ggplot(aes(P, col=targets)) + geom_freqpoly(binwidth=0.05)
res_RER_all %>% filter(tree == "red", direction=="down") %>% ggplot(aes(P, col=targets)) + geom_freqpoly(binwidth=0.05)

#SUBSETS TO TEST EFFECT

#trying subsets
vl_ext_sub = drop.tip(vl_ext, c("taeGut", "ficAlb", "corBra", "serCan", "geoFor","aptOwe", "aptHaa", "aptRow", "casCas", "rheAme", "rhePen"))
vl_ext_sub$edge.length[vl_ext_sub$edge.length > 1] = 1

rat_ext_sub = drop.tip(ratite_ext, c("aptOwe", "aptHaa", "aptRow", "casCas", "rheAme", "rhePen","taeGut", "ficAlb", "corBra", "serCan", "geoFor")) 
rat_ext_sub$edge.length[rat_ext_sub$edge.length > 1] = 1

fl_ext_sub = drop.tip(flightless_ext, c("aptOwe", "aptHaa", "aptRow", "casCas", "rheAme", "rhePen","taeGut", "ficAlb", "corBra", "serCan", "geoFor")) 
fl_ext_sub$edge.length[fl_ext_sub$edge.length > 1] = 1

RERext_sub <- RERconverge::getAllResiduals(aaSpTreesExt, useSpecies=vl_ext_sub$tip.label, transform = "sqrt", weighted = T, scale = T)

#make paths -- 7 total
vl_ext_paths_sub <- tree2Paths(vl_ext_sub, aaSpTreesExt)
ratite_ext_paths_sub <- tree2Paths(rat_ext_sub, aaSpTreesExt)
flightless_ext_paths_sub <- tree2Paths(fl_ext_sub, aaSpTreesExt)

### END MAKE PHENOTYPE TREES ### 

## RUN RER - 7 runs ###

res_vl_ext_sub <- correlateWithBinaryPhenotype(RERext_sub, vl_ext_paths_sub, min.sp=10, min.pos=3, weighted = "auto") %>% as_tibble(rownames="gene")
res_ratite_ext_sub <- correlateWithBinaryPhenotype(RERext_sub, ratite_ext_paths_sub, min.sp=10, min.pos=3, weighted = "auto") %>% as_tibble(rownames="gene")
res_flightless_ext_sub <- correlateWithBinaryPhenotype(RERext_sub, flightless_ext_paths_sub, min.sp=10, min.pos=3, weighted = "auto") %>% as_tibble(rownames="gene")

res_RER_all_sub <- bind_rows(list(vl_extsub = res_vl_ext_sub, rat_extsub = res_ratite_ext_sub, fl_extsub = res_flightless_ext_sub), .id = "run") %>% separate(run, into=c("targets", "tree")) %>% mutate(direction = ifelse(Rho < 0, "down", "up"))
write_tsv(res_RER_all_sub, path="all_RER_results_subsets.tsv")

res_RER_all_sub %>% mutate(direction = case_when(P < 0.01 & Rho < 0 ~ "down", P < 0.01 & Rho > 0 ~ "up", TRUE ~ "ns")) %>%
  group_by(targets, tree) %>% count(direction) %>% spread(direction, n) %>% mutate(prop_sig = (down + up) / (down + up + ns))



## OLD CODE ##

#test a random tree using foreground2Tree
targets<-c("anoDid", "nanHar", "strCam", "rheAme", "rhePen", "aptHaa", "aptOwe", "aptRow", "casCas", "droNov")
sample_tree <- vl_ext_sub %>% compute.brlen(0)
nontarget<-sample_tree$tip.label[!(sample_tree$tip.label %in% targets)]
testRandom <- foreground2Tree(sample(nontarget)[1:3], aaSpTreesExt, clade="weighted")
corRandom <- correlateWithBinaryPhenotype(RERext_sub, tree2Paths(testRandom,aaSpTreesExt))
corRandom %>% mutate(direction = case_when(P < 0.01 & Rho < 0 ~ "down", P < 0.01 & Rho > 0 ~ "up", TRUE ~ "ns")) %>% 
  count(direction) %>% spread(direction, n) %>% mutate(prop_sig = (down + up) / (down + up + ns))
