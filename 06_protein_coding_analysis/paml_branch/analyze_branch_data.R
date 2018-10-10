library(tidyverse)
library(devtools)
library(foreach)
library(doParallel)

#using my fork of RERconverge that fixes a few small bugs
if (!is.element("RERconverge", installed.packages()[,1])) {
  install_github("tsackton/RERconverge")
}

library(RERconverge)

## FUNCTIONS ##

make_trait_tree <- function(master.tree, target.tips, sp.to.drop, collapseAnc = FALSE) {
  #copied in part from RERconverge foreground2Tree function, but modifying tree input
  
  #check to make sure that nothing from target.tips is in sp.to.drop
  if (!is_empty(intersect(target.tips,sp.to.drop))) {
    stop("Target tips cannot overlap with species to drop")
  }
  
  #get list of tips to remove
  sp.to.remove = intersect(sp.to.drop, master.tree$tip.label)
  
  #get subtree by removing tips from master.tree corresponding to species to drop
  #also set all branch lengths to 0
  sub.tree <- drop.tip(master.tree, sp.to.remove) %>% compute.brlen(0)
  
  #get edges corresponding to tips
  if (!collapseAnc) {
    sub.tree$edge.length[RERconverge:::nameEdges(sub.tree) %in% target.tips] = 1
    names(sub.tree$edge.length) = RERconverge:::nameEdges(sub.tree)
  }
  else {
    tip.vals = rep(0, length(sub.tree$tip.label))
    names(tip.vals) = sub.tree$tip.label
    tip.vals[target.tips] = 1
    tmp = cbind(as.character(tip.vals))
    rownames(tmp) = names(tip.vals)
    tip.vals = tmp
    ancres = ancestral.pars(sub.tree, df <- as.phyDat(tip.vals, 
                                              type = "USER", levels = unique(as.character(tip.vals))), 
                                              type = "ACCTRAN")
    ancres = unlist(lapply(ancres, function(x) { x[2] }))
    internalVals = ancres
    evals = matrix(nrow = nrow(sub.tree$edge), ncol = 2)
    eres = ancres
    evals[, 1] = eres[sub.tree$edge[, 1]]
    evals[, 2] = eres[sub.tree$edge[, 2]]
    sub.tree$edge.length = evals[, 2] - evals[, 1]
    sub.tree$edge.length[sub.tree$edge.length < 1] = 0
  }
  
  #add weighted edges for internal nodes in clades
  edgeIndex = which(sub.tree$edge.length > 0)
  for (i in edgeIndex) {
    clade.edges = RERconverge:::getAllCladeEdges(sub.tree, i)
    clade.edges = unique(c(i, clade.edges))
    sub.tree$edge.length[clade.edges] = 1/length(clade.edges)
  }
  return(sub.tree)
}

make_RER_paths <- function(treeObj, trait.tree) {
  tree2Paths(trait.tree, treeObj)
}

count_sister_taxa <- function(tree, tips) {
  #set up pairs
  pairs<-combn(tips, 2, simplify = FALSE)
  lapply(pairs, is.monophyletic, phy=phy) %>% unlist %>% sum
}

get_random_targets <- function(number, tips, tree) {
  num_sister = 10
  while (num_sister > 0) {
    targets<-sample(tips, number)
    num_sister = count_sister_taxa(tree, targets)  
  }
  return(targets)
}

compute_results <- function(target_species, remove_species, tree_to_use) {
  trait_tree <- make_trait_tree(trees[[tree_to_use]], target_species, remove_species)
  tree_object <- get_tree_object(tree_to_use)
  rer_obj <- RERconverge::getAllResiduals(tree_object, cutoff=8e-06, useSpecies=trait_tree$tip.label, transform = "sqrt", weighted = T, scale = F)
  rer_res <- correlateWithBinaryPhenotype(rer_obj, make_RER_paths(tree_object, trait_tree), min.sp=10, min.pos=3, weighted = "auto") %>%
    as_tibble(rownames="gene") %>% 
    mutate(tree_version = tree_to_use) %>% 
    mutate(target_sp = paste0(sort(target_species), collapse="-"))
}

## END FUNCTIONS ##

## PRELIMINARY ##

#get hog to gene id
if (file.exists("hog_to_gene.Robj")) {
  load("hog_to_gene.Robj")
} else {
  hog_to_gene <- read_tsv("../../03_homology/HOG_final_alignment_seqids") %>% 
    mutate(hogid = as.integer(sub("HOG2_", "", HOG, fixed=T))) %>% 
    mutate(hogname = sub("HOG2_", "HOG", HOG, fixed = T)) %>%
    filter(Taxon == "galGal") %>%
    dplyr::select(-HOG, -Fasta_seqid, -Taxon, -Transcript) %>%
    group_by(hogid) %>%
    mutate(gene = paste0(Gene, sep="", collapse=";")) %>%
    dplyr::select(-Gene)
  save(hog_to_gene, file="hog_to_gene.Robj")
}

#get trees

extended_tree <- read.tree("../final_tree_ext_proteins.nwk") %>% 
  drop.tip(c("allMis", "chrPic", "anoCar", "egrGar"))
original_tree <- read.tree("../final_tree_ext_proteins.nwk") %>% 
  drop.tip(c("allMis", "chrPic", "anoCar", "egrGar", "nanAur", "nanBra", "nanHar", "uriPel"))
reduced_tree <- read.tree("../final_tree_ext_proteins.nwk") %>% 
  drop.tip(c("allMis", "chrPic", "anoCar", "egrGar", "nanAur", "nanBra", "nanHar", "uriPel", "anoDid"))

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

trees <- c("extended" = extended_tree, "original" = original_tree, "reduced" = reduced_tree)

#simple function to get tree object 

get_tree_object <- function(x) {
  #NO ERROR CHECKING
  if (x == "reduced") {
    aaSpTreesRed
  } else {
    aaSpTreesExt
  }
}

## END PRELIMINARY ##

##MAIN CODE TO RUN ANALYSIS AND WRITE OUT

## VARIABLE DEFINITIONS ##

target_species <- c("anoDid", "strCam", "rhePen")
remove_species <- c("rheAme", "droNov", "casCas", "aptHaa", "aptOwe", "aptRow")
tree_to_use <- "extended"

# RUN IN PARALLEL #

registerDoParallel(cores=4)

test_res<-foreach(ts=c("rhePen", "rheAme", "aptHaa", "aptRow"), .packages=c("tidyverse", "RERconverge")) %dopar%
  compute_results(c(ts, "anoDid", "strCam"), setdiff(c("rhePen", "rheAme", "aptHaa", "aptRow", "aptOwe", "droNov", "casCas"),ts),"extended")
  

