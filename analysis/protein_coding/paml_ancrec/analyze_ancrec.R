setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/paml_ancrec/")
library(data.table)
library(ape)
library(tidyverse)

#read in files
ancrec.parsed<-fread("gunzip -c paml_M0_parsed.txt.gz")
ancrec.treekey<-ancrec.parsed[,c("hog", "treenum", "species_tree"), with=FALSE]
hog_info<-read.table("../all_hog_info.tsv", sep="\t", header=T)

hog_info$has_species_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree]
hog_info$has_gene_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree == F]

ancrec.muts<-fread("gunzip -c paml_M0_parsedmuts.txt.gz")
names(ancrec.muts)<-c("hog", "tree", "branchnames", "branchpair", "ratite", "clean", "type", "mutclass")

ancrec.muts.merge<-merge(ancrec.muts, ancrec.treekey, by.x=c("hog", "tree"), by.y=c("hog", "treenum"), all=T)

#get tree
fulltree<-read.tree("../final_neut_tree.nwk")
tipdist<-cophenetic.phylo(fulltree)

#recaculate missing after removing species not in tree
tips_to_keep = fulltree$tip.label[fulltree$tip.label %in% names(hog_info)]
hog_info$missing_ct_clean = apply(hog_info[,tips_to_keep], 1, function(x) sum(is.na(x)))
hog_info$dup_ct_clean = apply(hog_info[,tips_to_keep], 1, function(x) sum(x[!is.na(x)]>1))

#merge
ancrec.muts.merge<-as.data.table(merge(ancrec.muts.merge, hog_info, by="hog"))

#ANALYSIS STARTS HERE##

#FUNCTIONS
prep_data <- function (DF, tipdist) {
  #this is a nasty function that includes everything!
  ancrec.1<-DF
  ancrec.1.sum<-as.data.frame.matrix(table(ancrec.1$branchpair, ancrec.1$mutclass))
  ancrec.1.key<-unique(ancrec.1[,c("branchpair", "branchnames", "ratite"), with=FALSE], by=c("branchpair", "branchnames", "ratite"))
  
  #add between clade and within clade info to key: neo-neo, paleo-paleo, neo-paleo
  palaeo<-c("rheAme", "rhePen", "droNov", "casCas", "aptRow", "aptOwe", "aptHaa", "strCam", "tinGut", "cryCin", "notPer", "eudEle")
  ancrec.1.key$branch1 = sub("-\\w+", "", ancrec.1.key$branchnames)
  ancrec.1.key$branch2 = sub("\\w+-", "", ancrec.1.key$branchnames)
  ancrec.1.key$dist = apply(ancrec.1.key, 1, function(x) if(x[4] %in% rownames(tipdist) & x[5] %in% rownames(tipdist)) { tipdist[x[4], x[5]] } else { return(NA) })
  
  ancrec.1.key$clade = "neo-palaeo"
  ancrec.1.key$clade[ancrec.1.key$branch1 %in% palaeo & ancrec.1.key$branch2 %in% palaeo] = "palaeo-palaeo"
  ancrec.1.key$clade[!(ancrec.1.key$branch1 %in% palaeo) & !(ancrec.1.key$branch2 %in% palaeo)] = "neo-neo"
  
  ancrec.1.sum<-merge(ancrec.1.sum, ancrec.1.key, by.x="row.names", by.y="branchpair")
  ancrec.1.sum$color="black"
  ancrec.1.sum$color[ancrec.1.sum$ratite==2]="red"
  ancrec.1.sum$broadconv = ancrec.1.sum$parallel+ancrec.1.sum$convergent
  ancrec.1.sum$broaddiv = ancrec.1.sum$divergent + ancrec.1.sum$regular
  ancrec.1.sum$broadratio = ancrec.1.sum$broadconv/(ancrec.1.sum$broaddiv+1)
  ancrec.1.sum$totalsubs = ancrec.1.sum$broadconv + ancrec.1.sum$broaddiv
  return(ancrec.1.sum)
}
run_lm_exp <- function(DF) {
  #linear model
  mod.full<-lm(broadconv ~ broaddiv*exp(dist) + as.factor(ratite==2), data=DF)
  mod.noratite<-lm(broadconv ~ broaddiv*exp(dist), data=DF)
  print(anova(mod.noratite,mod.full,test="LRT"))
  print(summary(mod.full))
}
run_lm <- function(DF) {
  #linear model
  mod.full<-lm(broadconv ~ broaddiv*dist + as.factor(ratite==2), data=DF)
  mod.noratite<-lm(broadconv ~ broaddiv*dist, data=DF)
  print(anova(mod.noratite,mod.full,test="LRT"))
  print(summary(mod.full))
}
run_lm_log <- function(DF) {
  #linear model
  mod.full<-lm(broadconv ~ broaddiv*log(dist) + as.factor(ratite==2), data=DF)
  mod.noratite<-lm(broadconv ~ broaddiv*log(dist), data=DF)
  print(anova(mod.noratite,mod.full,test="LRT"))
  print(summary(mod.full))
}

#not used
make_fig <- function(DF, file, distcut=0.15, pval=NA) {
  ancrec.1.sum <- DF
  pdf(file = file)
  plot(broadconv ~ broaddiv, col=color, pch=15, data=ancrec.1.sum[ancrec.1.sum$dist < distcut,], log="xy", xlab="# of Divergent Subs", ylab="# of Convergent Subs", bty="l")
  nonratite.mod<-lm(broadconv ~ 0 + broaddiv, data=ancrec.1.sum[ancrec.1.sum$ratite < 2 & ancrec.1.sum$dist < distcut,])
  plotpoints=data.frame(x=seq(1,10000,10), y=seq(1,10000,10)*coef(nonratite.mod))
  points(plotpoints, type='l')
  legend(x="topleft", legend=c("Ratite-ratite branch pair", "Other branch pairs"), pch=15, col=c("red", "black"), bty="n")
  if(!is.na(pval)) {
    legend("bottomright", legend=paste0("Ratite ANOVA P-value: ", pval), bty="n")
  }
  dev.off()
}
get_conv_num<-function(DF) {
  ancrec.conv<-subset(DF, (mutclass=="convergent" | mutclass=="parallel"))
  ancrec.all<-DF
  ancrec.all.gene<-as.data.frame.matrix(table(ancrec.all$hog, ancrec.all$ratite==2))
  names(ancrec.all.gene)=c("non", "rat")
  ancrec.all.gene$prop.ratite = ancrec.all.gene$rat/(ancrec.all.gene$non+ancrec.all.gene$rat)
  ancrec.conv.gene<-as.data.frame.matrix(table(ancrec.conv$hog, ancrec.conv$ratite==2))
  names(ancrec.conv.gene)=c("non", "rat")
  ancrec.test<-merge(ancrec.all.gene, ancrec.conv.gene, by="row.names")
  #conv is y
  ancrec.test2<-subset(ancrec.test, prop.ratite > 0)
  ancrec.test2$pval<-apply(ancrec.test2, 1, function(x) prop.test(as.numeric(x[6]),as.numeric(x[5])+as.numeric(x[6]),as.numeric(x[4]),alt="g")$p.value)
  ancrec.test2$qval <- p.adjust(ancrec.test2$pval, method="fdr")
  
  return(ancrec.test2)
}

#run on subsets

ancrec.1<-subset(ancrec.muts.merge, species_tree == TRUE & dup_ct_clean == 0 & clean=="clean" & type=="tip-tip" & missing_ct_clean<3, select=c("hog", "branchnames", "branchpair", "ratite", "mutclass"))
length(unique(ancrec.1$hog))
#prep data
ancrec.1.sum <- prep_data(ancrec.1, tipdist)
#get lm
run_lm(ancrec.1.sum)
run_lm_log(ancrec.1.sum)
run_lm_exp(ancrec.1.sum)

#gene tree, same filtering otherwise
ancrec.2<-subset(ancrec.muts.merge, species_tree == FALSE & dup_ct_clean == 0 & clean=="clean" & type=="tip-tip" & missing_ct_clean<3, select=c("hog", "branchnames", "branchpair", "ratite", "mutclass"))
length(unique(ancrec.2$hog))
#prep data
ancrec.2.sum <- prep_data(ancrec.2, tipdist)
#get lm
run_lm(ancrec.2.sum)
run_lm_log(ancrec.2.sum)
run_lm_exp(ancrec.2.sum)

#gene tree, strict filtering 
ancrec.3<-subset(ancrec.muts.merge, species_tree == FALSE & dup_ct_clean == 0 & clean=="clean" & type=="tip-tip" & missing_ct_clean == 0, select=c("hog", "branchnames", "branchpair", "ratite", "mutclass"))
length(unique(ancrec.3$hog))
#prep data
ancrec.3.sum <- prep_data(ancrec.3, tipdist)
#get lm
run_lm(ancrec.3.sum)
run_lm_log(ancrec.3.sum)
run_lm_exp(ancrec.3.sum)


##BELOW HERE NOT RUN
#make fig
#not run
#make_fig(ancrec.1.sum, "Figure2A.pdf", distcut=0.25, pval=round(unlist(run_lm(ancrec.1.sum)[6][2,]),3))
#specific genes
ratite.conv.genes<-get_conv_num(ancrec.1)
table(ratite.conv.genes$qval < 0.05)


#make fig
make_fig(ancrec.2.sum, "Figure2A-genetree.pdf", pval=round(unlist(run_lm(ancrec.2.sum)[6][2,]),3))
#specific genes
ratite.conv.genes.2<-get_conv_num(ancrec.2)
table(ratite.conv.genes.2$qval < 0.05)

#include internal, gene tree
ancrec.2<-subset(ancrec.muts.merge, species_tree == FALSE & dup_ct_clean == 0 & clean=="clean" & type=="tip-tip" & missing_ct_clean==0, select=c("hog", "branchnames", "branchpair", "ratite", "mutclass"))
length(unique(ancrec.2$hog))
#prep data
ancrec.2.sum <- prep_data(ancrec.2, tipdist)
#get lm
run_lm(ancrec.2.sum)
#make fig
make_fig(ancrec.2.sum, "Figure2A-genetree-strict.pdf", pval=round(unlist(run_lm(ancrec.2.sum)[6][2,]),3))
#specific genes
ratite.conv.genes.2<-get_conv_num(ancrec.2)
table(ratite.conv.genes.2$qval < 0.05)


#include internal, species tree
ancrec.2<-subset(ancrec.muts.merge, species_tree == TRUE & dup_ct_clean == 0 & clean=="clean" & type=="tip-tip" & missing_ct_clean==0, select=c("hog", "branchnames", "branchpair", "ratite", "mutclass"))
length(unique(ancrec.2$hog))
#prep data
ancrec.2.sum <- prep_data(ancrec.2, tipdist)
#get lm
run_lm(ancrec.2.sum)
#make fig
make_fig(ancrec.2.sum, "Figure2A-sptree-strict.pdf", pval=round(unlist(run_lm(ancrec.2.sum)[6][2,]),3))
#specific genes
ratite.conv.genes.2<-get_conv_num(ancrec.2)
table(ratite.conv.genes.2$qval < 0.05)
