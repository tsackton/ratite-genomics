setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/hyphy_bsrel/")
library(data.table)

#make hog<->galGal key
hogs<-read.table("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/homology/new_hog_list.txt")
ncbikey<-read.table("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/annotation/galGalAnnot/CGNC_Gallus_gallus_20161020.txt", header=F, sep="\t", comment.char="", col.names=c("cgnc", "ncbi", "ensembl", "sym", "descr", "sp"), quote="")
hogs.galgal<-subset(hogs, V4=="galGal")
hogs.galgal$hog=sub("HOG2_","",hogs.galgal$V1,fixed=T)
hogs.galgal<-merge(hogs.galgal, ncbikey, all.x=T, all.y=F, by.x="V3", by.y="ncbi")

#load tree key from paml_ancrec 
ancrec.parsed<-fread("gunzip -c ../paml_ancrec/ancrec_parsed.out.gz")
paml.treekey<-ancrec.parsed[,c("hog", "treenum", "species_tree"), with=FALSE]
paml.treekey$tree = paste0("tree", paml.treekey$treenum)

#load bs rel results
merged = read.table("bs_rel_parsed_full_2017-02-06.txt", fill=T)

#names
#rand1	tree1	10003	1	5	1	9	3	44	corBra:pygAde:colLiv:anaPla:Node45:eudEle	corBra:melUnd:falPer:picPub:pygAde:
names(merged)=c("class", "tree", "hog", "tsel.s", "nsel.s", "tsel.n", "nsel.n", "tnon", "nnon", "strict_branches", "nom_branches")

merged=merge(merged, paml.treekey, by.x=c("hog", "tree"), by.y=c("hog", "tree"), all=T)
merged$totbranch = merged$tnon + merged$nnon + merged$tsel.n + merged$nsel.n

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
write.table(unique(rerun$hog), col.names = F, row.names = F, quote=F, file="hogs_to_rerun_Feb2017")

merged$total_sel.s = merged$tsel.s + merged$nsel.s
merged$total_sel.n = merged$tsel.n +  merged$nsel.n
merged$target_prop.s = merged$tsel.s / merged$total_sel.s
table(merged$total_sel.s, merged$class)
table(merged$tsel.s, merged$class)

#rename
bsrel<-merged

bsrel$target_lin = bsrel$tsel.n + bsrel$tnon

bsrel.sp = subset(bsrel, target_lin > 0 & species_tree==TRUE & totbranch != 0 & !is.na(totbranch))
bsrel.gt = subset(bsrel, target_lin > 0 & species_tree==FALSE & totbranch != 0 & !is.na(totbranch))

table(bsrel.gt$target_prop.s == 1, bsrel.gt$class)
table(bsrel.sp$target_prop.s == 1, bsrel.sp$class)

barplot(table(bsrel.gt$tsel.s[bsrel.gt$nsel.s == 0 & bsrel.gt$tsel.s > 1], bsrel.gt$class[bsrel.gt$nsel.s == 0 & bsrel.gt$tsel.s > 1]), ylab="Number of Genes Uniquely Selected in >1 Target Lineages", legend=T)
barplot(table(bsrel.sp$tsel.s[bsrel.sp$nsel.s == 0 & bsrel.sp$tsel.s > 1], bsrel.sp$class[bsrel.sp$nsel.s == 0 & bsrel.sp$tsel.s > 1]), ylab="Number of Genes Uniquely Selected in >1 Target Lineages", legend=T)

#ratite specific selection
bsrel.ratite<-subset(bsrel.sp, class=="ratites")
bsrel.ratite$hog[bsrel.ratite$tsel.s>0 & bsrel.ratite$nsel.n==0]
hogs.galgal$sym[hogs.galgal$hog %in% bsrel.ratite$hog[bsrel.ratite$tsel.s>0 & bsrel.ratite$nsel.s==0]]

#vl specific selection
bsrel.vl<-subset(bsrel.sp, class=="vl")
bsrel.vl[bsrel.vl$tsel.s>0 & bsrel.vl$nsel.s==0,]

hogs.galgal$sym[hogs.galgal$hog %in% bsrel.vl$hog[bsrel.vl$tsel.s>0 & bsrel.vl$nsel.s==0]]

#thinking about ways to test for excess selection in ratites / vl

bsrel.ratite$prob_sel = (bsrel.ratite$total_sel.s / bsrel.ratite$totbranch)

compute_convergence <- function(prob = prob, target = target, total = total) {
  target<-rbinom(1, target, prob)
  nontarget<-rbinom(1, total-target, prob)
  if (nontarget == 0) {
    return(target)
  }
  else {
    return(0)
  }
}

exp_multiratite<-numeric(1000)

for (i in 1:1000) {
  res<-mapply(compute_convergence, bsrel.ratite$prob_sel, bsrel.ratite$target_lin, bsrel.ratite$totbranch)
exp_multiratite[i]<-sum(res > 1)/sum(res > 0)
}

#ratite specific selection
write.table(hogs.galgal$V3[hogs.galgal$hog %in% bsrel.ratite$hog[bsrel.ratite$tsel.s > 0 & bsrel.ratite$nsel.s == 0]], file="ratite_specific_selection_ncbi", sep="\t", row.names = F, col.names = F, quote=F)
write.table(hogs.galgal$sym[hogs.galgal$hog %in% bsrel.ratite$hog[bsrel.ratite$tsel.s > 0 & bsrel.ratite$nsel.s == 0]], file="ratite_specific_selection_sym", sep="\t", row.names = F, col.names = F, quote=F)               
