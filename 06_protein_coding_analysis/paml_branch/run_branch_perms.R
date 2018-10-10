setwd("~/Projects/birds/ratite_compgen/ratite-genomics/06_protein_coding_analysis/paml_branch/")
library(data.table)
library(tidyverse)
library(ape)

#set cores and perms for whole run, can be individually modified later in code
PERMS <- 1

##SETUP##

#read tree
phy<-read.tree(file="../final_neut_tree.nwk")

#read in hog_info key
#read in ancrec parsing key for species tree info

ancrec.parsed<-fread("gunzip -c ../paml_ancrec/paml_M0_parsed.txt.gz")
ancrec.treekey<-ancrec.parsed[,c("hog", "treenum", "species_tree"), with=FALSE]
hog_info<-read.table("../../03_homology/all_hog_info.tsv", sep="\t", header=T)

hog_info$has_species_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree]
hog_info$has_gene_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree == F]

#load hog <-> chicken info
hog_to_gene <- read.table("../../03_homology/HOG_final_alignment_seqids", header=T, stringsAsFactors =F)
hog_to_gene <- hog_to_gene %>% tbl_df %>%
  mutate(hog = as.integer(sub("HOG2_", "", HOG, fixed=T))) %>% 
  filter(Taxon == "galGal") %>%
  dplyr::select(-HOG, -Fasta_seqid, -Taxon, -Transcript) %>%
  group_by(hog) %>%
  mutate(gene = paste0(Gene, sep="", collapse=";")) %>%
  dplyr::select(-Gene)

##ANALYSIS - FUNCTIONS##

#setup - load and clean data
prep_data <- function(file, ancrec.treekey, hog_info, check_missing = F) {
  #load data -- dn
  read_line = paste0("gunzip -c ", file)
  dn<-fread(read_line, header=F, sep=",")
  names(dn)<-c("hog", "tree", "parent.node", "desc.node", "branch.id", "dn", "target")
  
  #add hog_info
  dn$tree = sub("tree", "", dn$tree, fixed=T)
  dn$tree = as.integer(dn$tree)
  dn<-merge(dn, ancrec.treekey, by.x=c("hog", "tree"), by.y=c("hog", "treenum"), all.x=T, all.y=F)
  dn<-merge(dn, hog_info, by.x="hog", by.y="hog", all=T)
  #get missing runs
  if (check_missing) {
    check_for_missing<-unique(subset(dn, is.na(dn), select=c("hog", "has_species_tree", "has_gene_tree")))
    write.table(check_for_missing$hog, file="branch_hogs_torerun.txt", quote=F, sep="", row.names = F, col.names =  T)
    print(nrow(check_for_missing))
   }
  return(dn)
}

subset_clean_data <- function(DF, missing_cutoff = 2, dup_cutoff = 0, use_sptree = TRUE) {
  dn.clean = subset(DF, dup_ct <= dup_cutoff & missing_ct <= missing_cutoff  & species_tree == use_sptree, select=c("hog", "parent.node", "desc.node", "branch.id", "dn", "target"))
  return(dn.clean)
}

#function to do vector projection
proj_vect <- function(genevec, sptree) {
  as.matrix(genevec) - sptree %*% t(sptree) %*% as.matrix(genevec)
}

#actually compute normalized stat
normalize_branch_stat <- function(DF, filter=TRUE) {
  dn.clean<-as.data.table(DF)
  #make unit vector
  dn.clean[,dn.length.bygene:=sqrt(sum(dn^2)), by=list(hog)]
  dn.clean$dn.unit.bygene=dn.clean$dn/dn.clean$dn.length.bygene
  #make average tree (average of all branches a tree appears on)
  dn.clean[,dn.average.tree:=mean(dn.unit.bygene, na.rm=T), by=list(branch.id)]
  #convert to unit vector (this will be different for each species tree configuration)
  dn.clean[,dn.unit.sptree:=dn.average.tree/sqrt(sum(dn.average.tree^2)), by=.(hog)]
  dn.clean[,dn.norm := proj_vect(dn.unit.bygene, dn.unit.sptree), by=.(hog)]
  if (filter) {
    branch_freqs<-as.data.frame(table(dn.clean$branch.id))
    dn.clean<-dn.clean[dn.clean$branch.id %in% branch_freqs$Var1[branch_freqs$Freq >= 500],]
  }
  return(dn.clean)
}

compute_results <- function (x, groupvar) {
  if (inherits(try(ans<-wilcox.test(x ~ groupvar, conf.int=TRUE),silent=TRUE),"try-error"))
    return(list(pval=NA_real_, est=NA_real_))
  else
    return(list(pval=ans$p.value, est=ans$estimate))
}

simple_perm <- function(perm, DF, target) {
  target <- enquo(target)
  key<-DF %>% dplyr::select(desc.node, !!target) %>% distinct()
  key$perm = dplyr::select(key, !!target) %>% pull() %>% sample()
  dn.perm <- inner_join(key, dn.default) %>% as.data.table
  dn.perm[,compute_results(dn.norm, perm), by=hog] %>% dplyr::select(est, pval, hog)
}

count_sister_taxa <- function(tree, tips) {
  #set up pairs
  pairs<-combn(tips, 2, simplify = FALSE)
  lapply(pairs, is.monophyletic, phy=phy) %>% unlist %>% sum
}

get_targets <- function(number, tips, tree) {
  num_sister = 10
  while (num_sister > 0) {
    targets<-sample(tips, number)
    num_sister = count_sister_taxa(tree, targets)  
  }
  return(targets)
}

random_species <- function(perm, DF, count, tips, tree, exclude = FALSE) {
  tips_to_use <- sort(get_targets(count, tips, tree))
  key<-DF %>% dplyr::select(desc.node) %>% distinct()
  key$perm = key$desc.node %in% tips_to_use
  if (exclude) {
    tips_to_exclude = tips[!tips %in% tips_to_use]
    key <- filter(key, !desc.node %in% tips_to_exclude)
  }
  tip_string <- paste0(tips_to_use, collapse = "-")
  dn.perm <- inner_join(key, dn.default) %>% as.data.table
  dn.perm[,compute_results(dn.norm, perm), by=hog] %>% dplyr::select(est, pval, hog) %>% mutate(tips = tip_string)
}

## ANALYSIS STARTS HERE ###

## DATA PREP ##

dn<-prep_data(file="aa_extended_withcorm_flightless.parsed.csv.gz", ancrec.treekey = ancrec.treekey, hog_info = hog_info, check_missing = T)
dn.default<-subset_clean_data(dn)
dn.default$target[grepl("aptHaa-aptOwe-aptRow-casCas-droNov-rheAme-rhePen", dn.default$desc.node, fixed=T)]=FALSE
dn.default$target[grepl("aptHaa-aptOwe-aptRow-casCas-droNov", dn.default$desc.node, fixed=T)]=FALSE
dn.default %>% filter(target==TRUE) %>% count(desc.node)

dn.default<-normalize_branch_stat(dn.default)

#real results
dn.default[,compute_results(dn.norm, target), by=hog] %>% dplyr::select(est, pval, hog) %>% write_tsv("target.real")

## PERMS ##
#do permutations in parallel
num_perms <- PERMS

#set up lists

simple_ratite<-vector("list", num_perms)
simple_vl<-vector("list", num_perms)

species_random<-vector("list", num_perms)
species_ratite<-vector("list", num_perms)
species_vl<-vector("list", num_perms)

for (i in 1:num_perms) {
  #run everything and make lists
  simple_ratite[[i]]<-simple_perm(DF=dn.default, target=ratite)
  simple_vl[[i]]<-simple_perm(DF=dn.default, target=vl)
  
  species_random[[i]]<-random_species(DF=dn.default, count=3, tips=all_clades, tree=phy)
  species_ratite[[i]]<-random_species(DF=dn.default, count=3, tips=ratite_clades, tree=phy)
  species_vl[[i]]<-random_species(DF=dn.default, count=3, tips=vl_clades, tree=phy)
}

simple_ratite %>% bind_rows(.id="perm") %>% write_tsv("simple_ratite.perms")
simple_vl %>% bind_rows(.id="perm") %>% write_tsv("simple_vl.perms")
simple_wb %>% bind_rows(.id="perm") %>% write_tsv("simple_wb.perms")

species_random %>% bind_rows(.id="perm") %>% write_tsv("species_random.perms")
species_ratite %>% bind_rows(.id="perm") %>% write_tsv("species_ratite.perms")
species_vl %>% bind_rows(.id="perm") %>% write_tsv("species_vl.perms")
