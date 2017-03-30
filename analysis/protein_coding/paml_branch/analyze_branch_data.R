setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/paml_branch/")
library(data.table)
library(qvalue)
library(tidyr)

#read in hog_info key
#read in ancrec parsing key for species tree info

ancrec.parsed<-fread("gunzip -c ../paml_ancrec/paml_M0_parsed.txt.gz")
ancrec.treekey<-ancrec.parsed[,c("hog", "treenum", "species_tree"), with=FALSE]
hog_info<-read.table("../all_hog_info.tsv", sep="\t", header=T)

hog_info$has_species_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree]
hog_info$has_gene_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree == F]

#make hog<->galGal key
hogs<-read.table("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/homology/new_hog_list.txt")
ncbikey<-read.table("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/annotation/galGalAnnot/CGNC_Gallus_gallus_20161020.txt", header=F, sep="\t", comment.char="", col.names=c("cgnc", "ncbi", "ensembl", "sym", "descr", "sp"), quote="")
hogs.galgal<-subset(hogs, V4=="galGal")
hogs.galgal$hog=as.integer(sub("HOG2_","",hogs.galgal$V1,fixed=T))
hogs.galgal<-merge(hogs.galgal, ncbikey, all.x=T, all.y=F, by.x="V3", by.y="ncbi")

#load data -- dn
dn<-fread("gunzip -c dn_parsed.txt", header=F, sep=",")
names(dn)<-c("hog", "tree", "parent.node", "desc.node", "branch.id", "dn", "ratite", "bop", "wb", "vl", "rand1", "rand2")

#add hog_info
dn$tree = sub("tree", "", dn$tree, fixed=T)
dn$tree = as.integer(dn$tree)
dn<-merge(dn, ancrec.treekey, by.x=c("hog", "tree"), by.y=c("hog", "treenum"), all.x=T, all.y=F)
dn<-merge(dn, hog_info, by.x="hog", by.y="hog", all=T)

#get missing runs
check_for_missing<-unique(subset(dn, is.na(dn), select=c("hog", "has_species_tree", "has_gene_tree")))
write.table(check_for_missing$hog, file="branch_reruns", quote=F, sep="", row.names = F, col.names =  F)

#clean up dn dataset -- start with just species tree, no dups, limited missing
missing_cutoff = 2
dn.clean = subset(dn, dup_ct == 0 & missing_ct <= missing_cutoff  & species_tree, select=c("hog", "parent.node", "desc.node", "branch.id", "dn"))
length(unique(dn.clean$hog))

#define clades -- fixed (vocal learners, palaeognaths, tinamous, ratites)
all_clades<-names(hog_info)[2:40]
ratite_clades<-c("aptHaa", "aptOwe", "aptRow", "rheAme", "rhePen", "strCam", "casCas", "droNov")
vl_clades<-c('calAnn', 'corBra', 'serCan', 'geoFor', 'melUnd', 'pseHum', 'taeGut', 'ficAlb')
tin_clades<-c("tinGut", "cryCin", "notPer", "eudEle")
palaeo_clades<-c(ratite_clades, tin_clades)
non_ratite_clades<-all_clades[!(all_clades %in% ratite_clades)]
non_ratite_palaeo<-palaeo_clades[!(palaeo_clades) %in% ratite_clades]

dn.clean$ratite=sapply(strsplit(dn.clean$desc.node, "-"), function(x) sum(x %in% ratite_clades)/length(x) == 1)
dn.clean$vl=sapply(strsplit(dn.clean$desc.node, "-"), function(x) sum(x %in% vl_clades)/length(x) == 1)
dn.clean$nrpalaeo=sapply(strsplit(dn.clean$desc.node, "-"), function(x) sum(x %in% non_ratite_palaeo)/length(x) == 1)

#make unit vector
dn.clean[,dn.length.bygene:=sqrt(sum(dn^2)), by=list(hog)]
dn.clean$dn.unit.bygene=dn.clean$dn/dn.clean$dn.length.bygene
#make average tree (average of all branches a tree appears on)
dn.clean[,dn.average.tree:=mean(dn.unit.bygene, na.rm=T), by=list(branch.id)]
#convert to unit vector (this will be different for each species tree configuration)
dn.clean[,dn.unit.sptree:=dn.average.tree/sqrt(sum(dn.average.tree^2)), by=.(hog)]

#function to do vector projection
proj_vect <- function(genevec, sptree) {
  as.matrix(genevec) - sptree %*% t(sptree) %*% as.matrix(genevec)
}
dn.clean[,dn.norm1 := proj_vect(dn,dn.unit.sptree), by=.(hog)]
dn.clean[,dn.norm2 := proj_vect(dn.unit.bygene, dn.unit.sptree), by=.(hog)]

#checks
branch_freqs<-as.data.frame(table(dn.clean$branch.id))

#exclude branches in less than 500 hogs
dn.clean<-dn.clean[dn.clean$branch.id %in% branch_freqs$Var1[branch_freqs$Freq >= 500],]
length(unique(dn.clean$hog))

#more checks
hist(dn.clean$dn.norm1, breaks=100)
hist(dn.clean$dn.norm2, breaks=100)

#normalization checks
high_freq_branches = branch_freqs$Var1[branch_freqs$Freq > 5000]
branch_key = unique(dn.clean[,c("branch.id","ratite","vl"), with=F])
plot_branch_dists = data.frame(dn.norm2=dn.clean$dn.norm2[dn.clean$branch.id %in% high_freq_branches], branch.id=dn.clean$branch.id[dn.clean$branch.id %in% high_freq_branches])
plot_branch_dists = merge(plot_branch_dists, branch_key, by.x="branch.id", by.y="branch.id")
plot_branch_dists$color = "black"
plot_branch_dists$color[plot_branch_dists$ratite]="red"
plot_branch_dists$color[plot_branch_dists$vl]="blue"
plot_branch_dists_colors = unique(plot_branch_dists[,c("branch.id", "color")])
boxplot(plot_branch_dists$dn.norm2 ~ plot_branch_dists$branch.id, outline=F, col=plot_branch_dists_colors$color)

#compute p-values
dn.clean[,ratite.p:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm1 ~ ratite)$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog)]


dn.clean[,vl.p:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm1 ~ vl)$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog)]

dn.pval<-unique(dn.clean[,c("hog", "ratite.p", "vl.p"), with=F])
length(dn.pval$hog)
summary(qvalue(dn.pval$ratite.p))
summary(qvalue(dn.pval$vl.p))

#plotting
#hog_to_plot = 1086
#with(dn.clean[dn.clean$hog==hog_to_plot,], plot(sort(dn.norm), col=ifelse(ratite[order(dn.norm)],"red", ifelse(nrpalaeo[order(dn.norm)], "black", ifelse(vl[order(dn.norm)], "blue", "gray50"))), pch=16))

#permutations
dn.perm=dn.clean[,.(hog,desc.node,dn.norm)]
dn.perm.desc_nodes=data.frame(desc.node=unique(dn.perm$desc.node),test_clade=F, stringsAsFactors = F)
nreps=1000
for (i in 1:nreps) {
  clades_okay = F
  while (!clades_okay) {
    test_clades = sample(all_clades, 8)
    clades_okay = sum(test_clades %in% ratite_clades) < 2 & sum(test_clades %in% vl_clades) < 2 
  }
  dn.perm.desc_nodes$test_clade=sapply(strsplit(dn.perm.desc_nodes$desc.node, "-"), function(x) sum(x %in%test_clades)/length(x) == 1)
 dn.perm$test_clade = dn.perm$desc.node %in% dn.perm.desc_nodes$desc.node[dn.perm.desc_nodes$test_clade]
 newcol=paste0("rep",i)
 dn.perm[,(newcol):= { if (inherits(try(ans<-wilcox.test(dn.norm ~ test_clade)$p.value,silent=TRUE),"try-error"))
   NA_real_
   else
     ans
 }, by=list(hog)]
}

dn.perm.pval=unique(dn.perm[,c(1,5:length(dn.perm)), with=F])
fwrite(dn.perm.pval, file="raw_dn_perm_out.csv")

#testing out plotting

#option1

plot(density(apply(dn.perm.pval, 2, function(x) sum(x < 0.05)/length(x))[2:length(dn.perm.pval)]), col="black", xlim=c(0,0.25), xlab="Fraction of tests with P < 0.05", las=1, bty="n", main="")
ratite.p.frac=sum(dn.pval$ratite.p[!is.na(dn.pval$ratite.p)] < 0.05)/sum(!is.na(dn.pval$ratite.p))
vl.p.frac=sum(dn.pval$vl.p[!is.na(dn.pval$vl.p)] < 0.05)/sum(!is.na(dn.pval$vl.p))
arrows(x0=vl.p.frac,y0=2,x1=vl.p.frac,y1=0, col="firebrick", lwd=2)
arrows(x0=ratite.p.frac,y0=2,x1=ratite.p.frac,y1=0, col="blue", lwd=2)
text(x=vl.p.frac, y=2.5, labels=c("Vocal Learners"))
text(x=ratite.p.frac, y=2.5, labels=c("Ratites"))

#option2
permdist<-cut(as.data.frame(dn.perm.pval)[,2], breaks=seq(0,1,0.01),labels=F)
plot(table(permdist), type="l", col=rgb(100,100,100,alpha=50,maxColorValue=255), ylim=c(0,700), xaxt="n", ylab="Count", las=1, xlab="P-value")
for (i in 3:length(dn.perm.pval)) {
  permdist<-cut(as.data.frame(dn.perm.pval)[,i], breaks=seq(0,1,0.01))
  points(table(permdist), type="l", col=rgb(100,100,100,alpha=50,maxColorValue=255))
  
}

ratite<-cut(dn.pval$ratite.p, breaks=seq(0,1,0.01), labels=F)
vl<-cut(dn.pval$vl.p, breaks=seq(0,1,0.01), labels=F)

lines(table(ratite), type="l", col="blue", lwd=3, lty="dashed")
lines(table(vl), type="l", col="firebrick", lwd=3, lty="dashed")
axis(1, labels=seq(0,1,0.2), at=seq(0,100,20))
legend("topright", legend=c("random", "ratite", "vocal learners"), col=c("gray50", "blue", "firebrick"), lwd=3, lty="dashed")

#option3
rand_sig_counts<-apply(dn.perm.pval[,2:length(dn.perm.pval)], 1, function(x) sum(p.adjust(x, method="fdr")<0.05))
plot(y=as.numeric(table(rand_sig_counts)), x=as.numeric(names(table(rand_sig_counts))), type="l", lwd=4, bty="n", xlim=c(0,400), xlab="# of Significant Genes", ylab="Freq", las=1)

vl.qvalue<-p.adjust(dn.pval$vl.p, method="fdr")
ratite.qvalue<-p.adjust(dn.pval$ratite.p, method="fdr")
vl.count<-sum(vl.qvalue[!is.na(vl.qvalue)]<0.05)
ratite.count<-sum(ratite.qvalue[!is.na(ratite.qvalue)]<0.05)
arrows(x0=vl.count,y0=500,x1=vl.count,y1=0, col="firebrick", lwd=2)
arrows(x0=ratite.count+25,y0=500,x1=ratite.count,y1=0, col="blue", lwd=2)
text(x=vl.count, y=600, labels=c("Vocal Learners"))
text(x=ratite.count+25, y=600, labels=c("Ratites"))

## VOCAL LEARNING GENES ###

vl.genes<-dn.pval[,.(hog,vl.p)]
vl.genes$vl.q<-p.adjust(vl.genes$vl.p, method="fdr")
vl.genes<-merge(vl.genes, hogs.galgal, by="hog", all.x=T, all.y=F)

#retain only those without NA in NCBI gene ID
vl.genes=subset(vl.genes, !is.na(V3) & !is.na(vl.p) & !is.na(vl.q))
#note this is higher because of a few duplications, which must have arisen due to some chicken dups/split gene models getting dropped from PAML alignments
sum(vl.genes$vl.q < 0.05)
vl.target<-vl.genes$V3[vl.genes$vl.q < 0.05]
vl.background<-vl.genes$V3
write.table(vl.target, file="vocal_learning_target_list.txt", sep="", quote=F, col.names = F, row.names = F)
write.table(vl.background, file="vocal_learning_background_list.txt", sep="", quote=F, col.names = F, row.names = F)

#get cor for pos/neg analysis
dn.clean[,vl.cor:= { 
  if (inherits(try(ans<-cor(dn.norm,as.numeric(vl),method="sp"),silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog)]

vl.cor<-unique(dn.clean[,c("hog", "vl.cor"), with=F])
vl.genes<-merge(vl.genes, vl.cor, by="hog", all.x=T, all.y=F)

write.table(vl.genes$V3[vl.genes$vl.q < 0.05 & vl.genes$vl.cor > 0], file="vocal_learning_accel_list.txt", sep="", quote=F, col.names = F, row.names = F)
write.table(vl.genes$V3[vl.genes$vl.q < 0.05 & vl.genes$vl.cor < 0], file="vocal_learning_decel_list.txt", sep="", quote=F, col.names = F, row.names = F)

head(vl.genes[order(vl.genes$vl.cor, decreasing=T),], n=25)
