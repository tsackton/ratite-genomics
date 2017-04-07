setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/hyphy_bsrel/")
library(data.table)
library(tidyr)
library(dplyr)

#load tree key from paml_ancrec 
ancrec.parsed<-fread("gunzip -c ../paml_ancrec/paml_M0_parsed.txt.gz")
ancrec.treekey<-ancrec.parsed[,c("hog", "treenum", "species_tree"), with=FALSE]
ancrec.treekey$tree = paste0("tree", paml.treekey$treenum)
hog_info<-read.table("../all_hog_info.tsv", sep="\t", header=T)
hog_info$has_species_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree]
hog_info$has_gene_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree == F]

#load hog <-> chicken info
hog_to_gene <- read.table("../HOG_final_alignment_seqids", header=T, stringsAsFactors =F)
hog_to_gene <- hog_to_gene %>% tbl_df %>%
  mutate(hog = as.integer(sub("HOG2_", "", HOG, fixed=T))) %>% 
  filter(Taxon == "galGal") %>%
  select(-HOG, -Fasta_seqid, -Taxon, -Transcript) %>%
  group_by(hog) %>%
  mutate(gene = paste0(Gene, sep="", collapse=";")) %>%
  select(-Gene)

##LOAD DATA## 
merged = read.table("bsrel_res_parsed_ratites_2017-04-06.txt", fill=T)
names(merged)=c("class", "tree", "hog", "tsel.s", "nsel.s", "tsel.n", "nsel.n", "tnon", "nnon", "strict_branches", "nom_branches")
merged=merge(merged, ancrec.treekey, by.x=c("hog", "tree"), by.y=c("hog", "tree"), all=T)

#work with tibbles
merged <- tbl_df(merged)
merged <- hog_info %>% select(hog, dup_ct, missing_ct) %>% right_join(., merged) %>% tbl_df %>% 
  mutate(totbranch = tnon + nnon + tsel.n + nsel.n) %>% 
  mutate(total_sel.s = tsel.s + nsel.s, total_sel.n = tsel.n + nsel.n, target_prop.s = tsel.s/total_sel.s, target_prop.n = tsel.n / total_sel.n, target_lin = tsel.n + tnon) %>% 
  select(-class)

## THIS CODE CHECKS FOR MISSING HOGS TO RERUN
#get missing and set up reruns
check_for_missing = subset(merged, class=="ratites", select=c("hog", "tree", "treenum", "species_tree", "totbranch"))

table(check_for_missing$treenum, check_for_missing$tree, useNA="ifany")
check_for_missing=subset(check_for_missing, !is.na(treenum))
table(check_for_missing$totbranch==0, useNA='ifany')
#missing set to rerun
rerun<-subset(check_for_missing, totbranch==0 | is.na(totbranch))
#rerun rule: if tree1 is a rerun & the species tree, rerun tree1 and tree2
#if tree1 is a rerun and not the species tree, no tree2 to rerun
#if tree2 is a rerun and not the species tree, no need to rerun tree1
#however because of checks in script, should be able to just run straight up
write.table(unique(rerun$hog), col.names = F, row.names = F, quote=F, file="hogs_to_rerun_April2017")
#####

#ANALYSIS 
bsrel<-merged %>% filter(dup_ct == 0, species_tree == TRUE, totbranch > 0, missing_ct <= 5)
bsrel %>% distinct(hog) %>% summarize(count=n())

barplot(table(bsrel$tsel.s[bsrel$nsel.s == 0 & bsrel$tsel.s >= 1]), ylab="Number of Genes Uniquely Selected in >1 Target Lineages", legend=T)

ratite_genes <- bsrel %>% filter(tsel.s > 1, nsel.s == 0) %>% select(hog) %>% inner_join(., hog_to_gene) %>% select(gene)
bsrel %>% filter(tsel.s > 1, nsel.s == 0) %>% select(hog) %>% inner_join(., hog_to_gene) %>% print.data.frame

all_genes <- bsrel %>% select(hog) %>% inner_join(., hog_to_gene) %>% select(gene)
write.table(ratite_genes, file="ratite_specific_selection_sp_relaxed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(all_genes, file="bsrel_hyphy_background_sp_relaxed", sep="\t", quote=F, row.names=F, col.names=F)

#thinking about ways to test for excess selection in ratites

bsrel$prob_sel = (bsrel$total_sel.s / bsrel$totbranch)

compute_convergence <- function(prob = prob, target = target, total = total) {
  target<-rbinom(1, target, prob)
  nontarget<-rbinom(1, total-target, prob)
  if (nontarget == 0) {
    return(target)
  }
  else {
    return(NA_real_)
  }
}

exp_multiratite<-list()

for (i in 1:1000) {
  res<-mapply(compute_convergence, bsrel$prob_sel, bsrel$target_lin, bsrel$totbranch)
exp_multiratite[[i]]<-as.data.frame(table(res), stringsAsFactors = F)
}

perm_res <- bind_rows(exp_multiratite, .id="rep") %>% spread(res, Freq, fill=0, drop=FALSE) %>% 
  rename(ct0 = `0`, ct1 = `1`, ct2 = `2`, ct3 = `3`, ct4=`4`, ct5=`5`) %>%
  mutate(selected = ct1+ct2+ct3+ct4+ct5, convergent = ct2+ct3+ct4+ct5, rel_conv = convergent/selected) %>% tbl_df


real_res <- as.data.frame(table(bsrel$tsel.s[bsrel$nsel.s==0])) %>% tbl_df %>% spread(Var1, Freq) %>%
  rename(ct0 = `0`, ct1 = `1`, ct2 = `2`, ct3 = `3`) %>%
  mutate(ct4=0, ct5=0, selected = ct1+ct2+ct3+ct4+ct5, convergent = ct2+ct3+ct4+ct5, rel_conv = convergent/selected) %>% tbl_df


1-(sum(real_res$selected > perm_res$selected)/1000)
1-(sum(real_res$convergent > perm_res$convergent)/1000)
1-(sum(real_res$rel_conv > perm_res$rel_conv)/1000)


bsrel %>% filter(tsel.s > 2, nsel.s == 0) %>% select(hog) %>% inner_join(., hog_to_gene) %>% arrange(gene) %>% print.data.frame
