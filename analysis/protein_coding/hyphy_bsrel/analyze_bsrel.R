setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/hyphy_bsrel/")
library(data.table)
ratites<-as.data.table(read.table("bsrel_res_parsed_ratites.txt", header=F, fill=T))
names(ratites)=c("class", "tree", "hog", "tsel", "nsel", "tnon", "nnon")
ratite.tree1<-subset(ratites, tree=="tree1")
ratite.tree1[,pval:={
  if (inherits(try(ans<-fisher.test(matrix(c(tsel,nsel,tnon,nnon),nrow=2))$p.value,silent=T),"try-error"))
    NA_real_
  else
    ans
  }, by=hog]



vl<-as.data.table(read.table("bsrel_res_parsed_vl.txt", header=F, fill=T))
names(vl)=c("class", "tree", "hog", "tsel", "nsel", "tnon", "nnon")
vl.tree1<-subset(vl, tree=="tree1")
vl.tree1[,pval:={
  if (inherits(try(ans<-fisher.test(matrix(c(tsel,nsel,tnon,nnon),nrow=2))$p.value,silent=T),"try-error"))
    NA_real_
  else
    ans
}, by=hog]

merged = read.table("merged_res.txt", fill=T)
names(merged)=c("class", "tree", "hog", "tsel", "nsel", "tnon", "nnon")
merged=subset(merged, tree=="tree1")
#looking at another approach
merged$total_sel = merged$tsel + merged$nsel
merged$target_prop = merged$tsel / merged$total_sel
plot(sort(merged$target_prop[merged$class=="ratites"]), type='l', col="red", ylab="Proportion of Selected Branches in Target Class")
lines(sort(merged$target_prop[merged$class=="vl"]), type="l", col="darkgreen")
lines(sort(merged$target_prop[merged$class=="rand1"]), type="l", col="gray")
lines(sort(merged$target_prop[merged$class=="rand2"]), type="l", col="gray")
lines(sort(merged$target_prop[merged$class=="wb"]), type="l", col="blue")
lines(sort(merged$target_prop[merged$class=="bop"]), type="l", col="purple")
legend("topleft", legend=c("'Random'", "Ratites", "Waterbirds", "Vocal Learners", "Birds of Prey"), col=c("gray", "red", "blue", "darkgreen", "purple"), lwd=2)

table(merged$tsel, merged$target_prop == 1, merged$class)

barplot(table(merged$tsel[merged$nsel == 0 & merged$tsel > 1], merged$class[merged$nsel == 0 & merged$tsel > 1]), ylab="Number of Genes Uniquely Selected in >1 Target Lineages")
legend("topleft", legend=c("2", "3", "4"), col=c("gray20", "gray50", "gray80"), pch=15, cex=2)

pos<-subset(merged, tsel>=3 & nsel==0)
pos<-merge(pos, hogs.galgal, by="hog", all.x=T, all.y=F)
