library(ape)

alltrees<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/paml/branch/omega_trees.out", header=F, sep="\t", stringsAsFactors=F)
names(alltrees)=c("hog", "tree")

alldata<-data.frame(taxa=character(0), omega=numeric(0), hog=numeric(0))

for (i in c(1:length(alltrees$hog))) {
  if (alltrees$tree[i] != "") {
    curtree = read.tree(text=alltrees$tree[i])
    curdata = data.frame(taxa=curtree$tip.label, omega=curtree$edge.length[1:length(curtree$tip.label)], hog=alltrees$hog[i])
    alldata<-rbind(alldata, curdata)    
  } 
}

write.table(alldata, file="tip.omega.allhogs", sep="\t", quote=F)
alldata$taxchar=apply(alldata[,1,drop=F], 1, nchar)

#remove non-species trees
alldata.clean=droplevels(subset(alldata, taxchar==6))

library(data.table)
alldata.clean<-as.data.table(alldata.clean)

alldata.long<-dcast(alldata.clean, hog ~ taxa, value.var="omega")
alldata.normalized<-as.data.frame(scale(alldata.long[,2:40, with=F], scale=F))
alldata.normalized$accel.p = apply(alldata.normalized, 1, function(x) wilcox.test(x[c(1,2,6,7,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,30,31,32,33,34,35,36,37,38,39)], x[c(3,4,5,8,14,27,28,29)], alternative="less")$p.value)
alldata.normalized$accel.p2 = apply(alldata.normalized, 1, function(x) wilcox.test(x[c(1,2,3,4,5,6,7,8,9,12,13,14,15,16,17,18,20,21,24,26,27,28,29,30,31,32,36,37,38,33)], x[c(39,34,35,10,22,23,11,19,25)], alternative="less")$p.value)

alldata.normalized$decel.p = apply(alldata.normalized, 1, function(x) wilcox.test(x[c(1,2,6,7,9,10,11,12,13,15,16,17,18,19,20,21,23,24,25,26,30,31,32,33,34,35,36,37,38,39)], x[c(3,4,5,8,14,27,28,29)], alternative="greater")$p.value)
alldata.normalized$est = apply(alldata.normalized, 1, function(x) wilcox.test(x[c(1,2,6,7,9,10,11,12,13,15,16,17,18,19,20,21,23,24,25,26,30,31,32,33,34,35,36,37,38,39)], x[c(3,4,5,8,14,27,28,29)],conf.int=TRUE)$estimate)

#really no evidence for rate shifts

#unnormalized
head(alldata.long)
alldata.long$accel.p = apply(alldata.long[,c(2:40), with=F], 1, function(x) wilcox.test(as.numeric(x[c(1,2,6,7,9,10,11,12,13,15,16,17,18,19,20,21,23,24,25,26,30,31,32,33,34,35,36,37,38,39)]), as.numeric(x[c(3,4,5,8,14,27,28,29)]), alternative="less")$p.value)


hist(alldata.normalized$accel.p, breaks=30, col="red", cex.axis=2, las=1, xlab="", ylab="", main="")


hist(alldata.normalized$decel.p, breaks=30, col="blue", cex.axis=2, las=1, xlab="", ylab="", main="")


plot(density(as.numeric(alldata.normalized[7,c(1,2,6,7,9,10,11,12,13,15,16,17,18,19,20,21,23,24,25,26,30,31,32,33,34,35,36,37,38,39)]), na.rm=T), xlim=c(-1,1), bty="n", main="", xlab="", ylab="", las=1, lwd=2, col="blue")
points(y=rep(0.2,8), x=alldata.normalized[7,c(3,4,5,8,14,27,28,29)], pch=16, col="red")
  


