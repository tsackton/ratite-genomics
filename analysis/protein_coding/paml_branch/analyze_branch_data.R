setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/paml_branch/")
library(data.table)
library(qvalue)

w<-fread("gunzip -c dn_parsed.txt", header=F, sep=",")
names(w)<-c("hog", "tree", "parent.node", "desc.node", "branch.id", "w", "ratite", "bop", "wb", "vl", "rand1", "rand2")

#subset to keep only branches present frequently
branchfreq<-as.data.frame(table(w$branch.id))
w.clean<-subset(w, branch.id %in% branchfreq$Var1[branchfreq$Freq>1])

#normalize by branch
w.clean[,w.sum.bygene:=sum(w), by=list(hog,tree)]
w.clean$w.norm.bygene=w.clean$w/w.clean$w.sum.bygene
#w.norm.bygene should add to to 1 for each hog/tree -- it does
#now scale to be relative rate for each branch
w.clean[,w.mean.bybranch:=mean(w.norm.bygene,na.rm=T), by=list(branch.id,tree)]
w.clean$w.norm=w.clean$w.norm.bygene - w.clean$w.mean.bybranch

w.clean[,ratite.pg:= { 
  if (inherits(try(ans<-wilcox.test(w.norm[ratite], w.norm[!ratite], alternative="greater")$p.value,silent=TRUE),"try-error"))
       NA_real_
  else
       ans
}, by=list(hog,tree)]

w.clean[,ratite.pl:= { 
  if (inherits(try(ans<-wilcox.test(w.norm[ratite], w.norm[!ratite], alternative="l")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,ratite.p:= { 
  if (inherits(try(ans<-wilcox.test(w.norm ~ ratite)$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,ratite.cor:= { 
  if (inherits(try(ans<-cor(w.norm,as.numeric(ratite),method="sp"),silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]


w.clean[,wb.pg:= { 
  if (inherits(try(ans<-wilcox.test(w.norm[wb], w.norm[!wb], alternative="greater")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,wb.pl:= { 
  if (inherits(try(ans<-wilcox.test(w.norm[wb], w.norm[!wb], alternative="l")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,wb.p:= { 
  if (inherits(try(ans<-wilcox.test(w.norm ~ wb)$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]


w.clean[,rand1.pg:= { 
  if (inherits(try(ans<-wilcox.test(w.norm[rand1], w.norm[!rand1], alternative="greater")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,rand1.pl:= { 
  if (inherits(try(ans<-wilcox.test(w.norm[rand1], w.norm[!rand1], alternative="l")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,rand1.p:= { 
  if (inherits(try(ans<-wilcox.test(w.norm ~ rand1)$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,rand2.p:= { 
  if (inherits(try(ans<-wilcox.test(w.norm ~ rand2)$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,vl.p:= { 
  if (inherits(try(ans<-wilcox.test(w.norm ~ vl)$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,bop.p:= { 
  if (inherits(try(ans<-wilcox.test(w.norm ~ bop)$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]


w.clean[,rand1.cor:= { 
  if (inherits(try(ans<-cor(w.norm,as.numeric(rand1),method="sp"),silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,bop.cor:= { 
  if (inherits(try(ans<-cor(w.norm,as.numeric(bop),method="sp"),silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,rand2.cor:= { 
  if (inherits(try(ans<-cor(w.norm,as.numeric(rand2),method="sp"),silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,vl.cor:= { 
  if (inherits(try(ans<-cor(w.norm,as.numeric(vl),method="sp"),silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,wb.cor:= { 
  if (inherits(try(ans<-cor(w.norm,as.numeric(wb),method="sp"),silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]





#crude hack to get tinamou branches
w.clean$tinamou=sapply(strsplit(w.clean$desc.node, "-"), function(x) sum(x %in% c("tinGut", "cryCin", "notPer", "eudEle"))/length(x))


w.clean[,tinamou.pg:= { 
  if (inherits(try(ans<-wilcox.test(w.norm[tinamou==1], w.norm[tinamou<1], alternative="greater")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,tinamou.pl:= { 
  if (inherits(try(ans<-wilcox.test(w.norm[tinamou==1], w.norm[tinamou<1], alternative="l")$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.clean[,tinamou.p:= { 
  if (inherits(try(ans<-wilcox.test(w.norm[tinamou==1], w.norm[tinamou<1])$p.value,silent=TRUE),"try-error"))
    NA_real_
  else
    ans
}, by=list(hog,tree)]

w.pval<-unique(w.clean[tree=="tree1",c("hog","tree","ratite.pl","ratite.pg", "ratite.p", "rand1.pl", "rand1.pg", "rand1.p", "wb.pl", "wb.pg", "wb.p", "tinamou.pl", "tinamou.pg", "tinamou.p", "rand2.p", "vl.p", "bop.p"), with=F])

w.cor<-unique(w.clean[tree=="tree1", c("hog", "ratite.cor", "rand1.cor", "bop.cor", "rand2.cor", "wb.cor", "vl.cor"), with=F])

#plotting
hist(w.clean$w.norm[w.clean$hog==1293 & !w.clean$ratite], breaks=20, xlim=c(-0.1,0.1))
points(w.clean$w.norm[w.clean$hog==1293 & w.clean$ratite], rep(1,sum(w.clean$ratite[w.clean$hog==1293])), col="red", pch=16, cex=1.1)

w.cor$bop.scale=scale(w.cor$bop.cor, center=T, scale=F)
w.cor$wb.scale=scale(w.cor$wb.cor, center=T, scale=F)
w.cor$vl.scale=scale(w.cor$vl.cor, center=T, scale=F)
w.cor$rand1.scale=scale(w.cor$rand1.cor, center=T, scale=F)
w.cor$rand2.scale=scale(w.cor$rand2.cor, center=T, scale=F)
w.cor$ratite.scale=scale(w.cor$ratite.cor, center=T, scale=F)

plot(density(w.cor$rand1.scale, na.rm=T), col="black", lwd=2, main="", xlab="Normalized Correlation Coefficient")
lines(density(w.cor$rand2.scale, na.rm=T), col="black", lwd=2)
lines(density(w.cor$ratite.scale, na.rm=T), col="blue", lwd=2, lty="dashed")
#lines(density(w.cor$bop.scale, na.rm=T), col="red", lwd=2, lty="dashed")
#lines(density(w.cor$wb.scale, na.rm=T), col="purple", lwd=2, lty="dashed")
lines(density(w.cor$vl.scale, na.rm=T), col="green", lwd=2, lty="dashed")
abline(v=0)
legend("topright", legend=c("random", "ratite", "vocal learners"), col=c("black", "blue", "green"), lwd=2, lty="dashed")

hist(qvalue(w.pval$ratite.p))
summary(qvalue(w.pval$ratite.p))
hist(qvalue(w.pval$rand1.p))
summary(qvalue(w.pval$rand1.p))
hist(qvalue(w.pval$rand2.p))
summary(qvalue(w.pval$rand2.p))

hist(qvalue(w.pval$bop.p))
hist(qvalue(w.pval$wb.p))

hist(qvalue(w.pval$vl.p))
summary(qvalue(w.pval$vl.p))

w.pval$ratite.padj=p.adjust(w.pval$ratite.p, method="fdr")
w.pval$vl.padj=p.adjust(w.pval$vl.p, method="fdr")

vl.hogs<-w.pval$hog[w.pval$vl.padj<0.05]
vl.acc<-w.cor$hog[w.cor$vl.cor>0]
vl.hogs.acc<-vl.hogs[vl.hogs %in% vl.acc]

vl.hogs.acc.galgal<-hogs.galgal[hogs.galgal$hog %in% vl.hogs.acc,]
write.table(vl.hogs.acc.galgal$V3, file="vl.ncbi.acc", quote=F, col.names=F, row.names=F)
write.table(hogs.galgal$V3, file="branch.background", quote=F, col.names=F, row.names=F)
