setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/paml_branch/")
library(data.table)
library(qvalue)

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

dn<-fread("gunzip -c dn_parsed.txt", header=F, sep=",")
names(dn)<-c("hog", "tree", "parent.node", "desc.node", "branch.id", "dn", "ratite", "bop", "wb", "vl", "rand1", "rand2")

#add species tree info to dn
dn<-merge(dn, paml.treekey, by=c("hog", "tree"))

#subset to keep only branches present frequently
branchfreq<-as.data.frame(table(dn$branch.id))
dn.clean<-subset(dn, branch.id %in% branchfreq$Var1[branchfreq$Freq>1])


#normalize by branch
dn.clean[,dn.sum.bygene:=sum(dn), by=list(hog,species_tree)]
dn.clean$dn.norm.bygene=dn.clean$dn/dn.clean$dn.sum.bygene
#dn.norm.bygene should add to to 1 for each hog/tree -- it does
#check
dn.clean[,dn.check:=sum(dn.norm.bygene), by=list(hog,species_tree)]
#now scale to be relative rate for each branch
dn.clean[,dn.mean.bybranch:=mean(dn.norm.bygene,na.rm=T), by=list(branch.id,species_tree)]
dn.clean$dn.norm=dn.clean$dn.norm.bygene - dn.clean$dn.mean.bybranch

dn.clean[,ratite.pg:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm[ratite], dn.norm[!ratite], alternative="greater")$p.value,silent=TRUE),"try-error"))
       NA_real_
  else
       ans
}, by=list(hog,species_tree)]

dn.clean[,ratite.pl:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm[ratite], dn.norm[!ratite], alternative="l")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,ratite.p:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm ~ ratite)$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,ratite.cor:= { 
  if (inherits(try(ans<-cor(dn.norm,as.numeric(ratite),method="sp"),silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,rand1.pg:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm[rand1], dn.norm[!rand1], alternative="greater")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,rand1.pl:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm[rand1], dn.norm[!rand1], alternative="l")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,rand1.p:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm ~ rand1)$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]


dn.clean[,rand1.cor:= { 
  if (inherits(try(ans<-cor(dn.norm,as.numeric(rand1),method="sp"),silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,rand2.pg:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm[rand2], dn.norm[!rand2], alternative="greater")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,rand2.pl:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm[rand2], dn.norm[!rand2], alternative="l")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,rand2.p:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm ~ rand2)$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,rand2.cor:= { 
  if (inherits(try(ans<-cor(dn.norm,as.numeric(rand2),method="sp"),silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,vl.pg:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm[vl], dn.norm[!vl], alternative="greater")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,vl.pl:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm[vl], dn.norm[!vl], alternative="l")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]


dn.clean[,vl.p:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm ~ vl)$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,vl.cor:= { 
  if (inherits(try(ans<-cor(dn.norm,as.numeric(vl),method="sp"),silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

#crude hack to get tinamou branches
dn.clean$tinamou=sapply(strsplit(dn.clean$desc.node, "-"), function(x) sum(x %in% c("tinGut", "cryCin", "notPer", "eudEle"))/length(x))

dn.clean[,tinamou.pg:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm[tinamou==1], dn.norm[tinamou<1], alternative="greater")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,tinamou.pl:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm[tinamou==1], dn.norm[tinamou<1], alternative="l")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,tinamou.p:= { 
  if (inherits(try(ans<-wilcox.test(dn.norm[tinamou==1], dn.norm[tinamou<1])$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.clean[,tinamou.cor:= { 
  if (inherits(try(ans<-cor(dn.norm,as.numeric(tinamou==1),method="sp"),silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,species_tree)]

dn.pval.st<-unique(dn.clean[species_tree==T,c("hog","species_tree","ratite.pl","ratite.pg", "ratite.p", "rand1.pl", "rand1.pg", "rand1.p", "tinamou.pl", "tinamou.pg", "tinamou.p", "rand2.pl", "rand2.pg", "rand2.p", "vl.pl", "vl.pg", "vl.p"), with=F], by=c("hog", "species_tree"))

dn.pval.gt<-unique(dn.clean[species_tree==F,c("hog","species_tree","ratite.pl","ratite.pg", "ratite.p", "rand1.pl", "rand1.pg", "rand1.p", "tinamou.pl", "tinamou.pg", "tinamou.p", "rand2.pl", "rand2.pg", "rand2.p", "vl.pl", "vl.pg", "vl.p"), with=F], by=c("hog", "species_tree"))

dn.cor.st<-unique(dn.clean[species_tree==T,c("hog", "species_tree", "ratite.cor", "rand1.cor", "rand2.cor", "vl.cor", "tinamou.cor"), with=F], by=c("hog", "species_tree"))

dn.cor.gt<-unique(dn.clean[species_tree==F, c("hog", "species_tree", "ratite.cor", "rand1.cor", "rand2.cor", "vl.cor", "tinamou.cor"), with=F], by=c("hog", "species_tree"))

#qvalue analysis - species tree

hist(qvalue(dn.pval.st$ratite.p))
summary(qvalue(dn.pval.st$ratite.p))
hist(qvalue(dn.pval.st$rand1.p))
summary(qvalue(dn.pval.st$rand1.p))
hist(qvalue(dn.pval.st$rand2.p))
summary(qvalue(dn.pval.st$rand2.p))
hist(qvalue(dn.pval.st$vl.p))
summary(qvalue(dn.pval.st$vl.p))

#final plot draft
ratite<-cut(dn.pval.st$ratite.p, breaks=seq(0,1,0.01), include.lowest=F, right=T, labels=F)
rand1<-cut(dn.pval.st$rand1.p, breaks=seq(0,1,0.01), include.lowest=F, right=T, labels=F)
rand2<-cut(dn.pval.st$rand2.p, breaks=seq(0,1,0.01), include.lowest=F, right=T, labels=F)
vl<-cut(dn.pval.st$vl.p, breaks=seq(0,1,0.01), include.lowest=F, right=T, labels=F)

plot(table(vl), type="l", col="firebrick", xlab="P-value", ylab="Count", bty="l", las=1, xaxt="n", lwd=3, lty="dashed")
lines(table(ratite), type="l", col="blue", lwd=3, lty="dashed")
lines(table(rand1), type="l", col="black", lwd=3, lty="dashed")
lines(table(rand2), type="l", col="black", lwd=3, lty="dashed")
axis(1, labels=seq(0,1,0.2), at=seq(0,100,20))
legend("topright", legend=c("random", "ratite", "vocal learners"), col=c("black", "blue", "firebrick"), lwd=3, lty="dashed")

#qvalue analysis - gene tree

hist(qvalue(dn.pval.gt$ratite.p))
summary(qvalue(dn.pval.gt$ratite.p))
hist(qvalue(dn.pval.gt$rand1.p))
summary(qvalue(dn.pval.gt$rand1.p))
hist(qvalue(dn.pval.gt$rand2.p))
summary(qvalue(dn.pval.gt$rand2.p))
hist(qvalue(dn.pval.gt$vl.p))
summary(qvalue(dn.pval.gt$vl.p))

#cor analysis - sp tree

dn.cor.st$vl.scale=scale(dn.cor.st$vl.cor, center=T, scale=F)
dn.cor.st$rand1.scale=scale(dn.cor.st$rand1.cor, center=T, scale=F)
dn.cor.st$rand2.scale=scale(dn.cor.st$rand2.cor, center=T, scale=F)
dn.cor.st$ratite.scale=scale(dn.cor.st$ratite.cor, center=T, scale=F)

plot(density(dn.cor.st$rand1.scale, na.rm=T), col="black", lwd=2, main="", xlab="Normalized Correlation Coefficient")
lines(density(dn.cor.st$rand2.scale, na.rm=T), col="black", lwd=2)
lines(density(dn.cor.st$ratite.scale, na.rm=T), col="blue", lwd=2, lty="dashed")
lines(density(dn.cor.st$vl.scale, na.rm=T), col="firebrick", lty="dashed", lwd=2)
abline(v=0)
legend("topright", legend=c("random", "ratite", "vocal learners"), col=c("black", "blue", "firebrick"), lwd=2, lty="dashed")

#cor analysis - gene tree

dn.cor.gt$vl.scale=scale(dn.cor.gt$vl.cor, center=T, scale=F)
dn.cor.gt$rand1.scale=scale(dn.cor.gt$rand1.cor, center=T, scale=F)
dn.cor.gt$rand2.scale=scale(dn.cor.gt$rand2.cor, center=T, scale=F)
dn.cor.gt$ratite.scale=scale(dn.cor.gt$ratite.cor, center=T, scale=F)

plot(density(dn.cor.gt$rand1.scale, na.rm=T), col="black", lwd=2, main="", xlab="Normalized Correlation Coefficient")
lines(density(dn.cor.gt$rand2.scale, na.rm=T), col="black", lwd=2)
lines(density(dn.cor.gt$ratite.scale, na.rm=T), col="blue", lwd=2, lty="dashed")
lines(density(dn.cor.gt$vl.scale, na.rm=T), col="firebrick", lty="dashed", lwd=2)
abline(v=0)
legend("topright", legend=c("random", "ratite", "vocal learners"), col=c("black", "blue", "firebrick"), lwd=2, lty="dashed")

#analyze vocal learning results
vl.hogs<-data.frame(hog=dn.pval.st[p.adjust(dn.pval.st$vl.p, method="fdr") < 0.05, hog],acc=F)
vl.hogs$acc[vl.hogs$hog %in% dn.cor.st$hog[dn.cor.st$vl.cor > 0]]=T
vl.hogs<-merge(vl.hogs, hogs.galgal, by.x="hog", by.y="hog")

#make lists of protein ids
write.table(unique(vl.hogs$V2[vl.hogs$acc==T]), file="vl.ncbi.acc", quote=F, col.names=F, row.names=F)
write.table(unique(vl.hogs$V2[vl.hogs$acc==F]), file="vl.ncbi.decel", quote=F, col.names=F, row.names=F)
write.table(unique(hogs.galgal$V2), file="vl.background", quote=F, col.names=F, row.names=F)

#compare to relax
source("../hyphy_relax/analyze_hyphy_relax.R")

comp<-merge(dn.pval.st, relax.tips.wide, by="hog")
fisher.test(table(comp$sig.vl != 0, factor(p.adjust(comp$vl.p, method="fdr") < 0.05, levels=c("FALSE", "TRUE"))))

comp2<-merge(dn.cor.st, relax.tips.wide, by="hog")
head(comp2)
boxplot(comp2$vl.cor ~ comp2$sig.vl, notch=T, outline=F)

table(comp$sig.rand != 0)
table(comp$sig.vl != 0)
table(comp$sig.ratite != 0)

boxplot(-log10(comp$ratite.p) ~ comp$sig.ratite != 0, outline=F, notch=T, ylim=c(0,5))
boxplot(-log10(comp$vl.pg) ~ comp$sig.vl, outline=F, notch=T, ylim=c(0,5))
boxplot(-log10(comp$vl.pl) ~ comp$sig.vl, outline=F, notch=T, ylim=c(0,5))

boxplot(-log10(comp$rand1.p) ~ comp$sig.rand != 0, outline=F, notch=T, ylim=c(0,5))
