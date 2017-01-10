setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/paml_ancrec/")
library(data.table)
ancrec<-fread("gunzip -c parsed_mutations.txt.gz")
names(ancrec)<-c("hog", "tree", "branchpair", "ratite", "clean", "type", "mutclass")

#make clean subset
ancrec.clean<-subset(ancrec, clean=="clean" & tree=="1")

#branchpair <-> ratite key
branchpair.key<-unique(ancrec.clean[,c("branchpair", "ratite"), with=FALSE])

#branchpair <-> hog count
branchpair.hog<-unique(ancrec.clean[,c("hog", "branchpair"), with=FALSE])
branchpair.count<-as.data.frame(table(branchpair.hog$branchpair))
branchpair.key<-merge(branchpair.key, branchpair.count, by.x="branchpair", by.y="Var1")

#all
ancrec.all.sum<-as.data.frame.matrix(table(ancrec.clean$branchpair, ancrec.clean$mutclass))
ancrec.all.sum$total = ancrec.all.sum$convergent + ancrec.all.sum$divergent + ancrec.all.sum$parallel + ancrec.all.sum$regular
ancrec.all.sum<-merge(ancrec.all.sum, branchpair.key, by.x="row.names", by.y="branchpair", all.x=T, all.y=F)
ancrec.all.sum<-subset(ancrec.all.sum, Freq > 10 & total > 100)
ancrec.all.sum$plotcol = "red"
ancrec.all.sum$plotcol[ancrec.all.sum$ratite==1] = "salmon"
ancrec.all.sum$plotcol[ancrec.all.sum$ratite==2] = "black"

plot(ancrec.all.sum$convergent ~ ancrec.all.sum$divergent, pch=1, cex=1, col=ancrec.all.sum$plotcol, log="xy")

ancrec.all.sum$ratio = log2(ancrec.all.sum$convergent / (ancrec.all.sum$divergent+1))
hist(ancrec.all.sum$ratio[ancrec.all.sum$ratite < 2], breaks=50, main="All Branches", xlab="Convergent/Divergent Ratio", xlim=c(-10,0))
points(x=ancrec.all.sum$ratio[ancrec.all.sum$ratite==2], y=rep(25, sum(ancrec.all.sum$ratite==2)), pch=16, col="red", cex=1.5)

#just tips
ancrec.clean.tips=subset(ancrec.clean, type=="tip-tip")
ancrec.tips.sum<-as.data.frame.matrix(table(ancrec.clean.tips$branchpair, ancrec.clean.tips$mutclass))
ancrec.tips.sum$total = ancrec.tips.sum$convergent + ancrec.tips.sum$divergent + ancrec.tips.sum$parallel + ancrec.tips.sum$regular
ancrec.tips.sum<-merge(ancrec.tips.sum, branchpair.key, by.x="row.names", by.y="branchpair", all.x=T, all.y=F)
ancrec.tips.sum<-subset(ancrec.tips.sum, Freq > 10 & total > 100)
ancrec.tips.sum$plotcol = "red"
ancrec.tips.sum$plotcol[ancrec.tips.sum$ratite==1] = "salmon"
ancrec.tips.sum$plotcol[ancrec.tips.sum$ratite==2] = "black"

#plot(ancrec.tips.sum$convergent ~ ancrec.tips.sum$divergent, pch=1, cex=1, col=ancrec.tips.sum$plotcol, log="xy")

ancrec.tips.sum$ratio = log2(ancrec.tips.sum$convergent / (ancrec.tips.sum$divergent+1))
hist(ancrec.tips.sum$ratio[ancrec.tips.sum$ratite < 2], breaks=50, main="Tips Only", xlab="Convergent/Divergent Ratio", xlim=c(-10,0))
points(x=ancrec.tips.sum$ratio[ancrec.tips.sum$ratite==2], y=rep(7, sum(ancrec.tips.sum$ratite==2)), pch=16, cex=1.5, col="red")

#just convergent sites
ancrec.conv<-subset(ancrec.clean, mutclass=="convergent" & type=="tip-tip")
ancrec.conv.gene<-as.data.frame.matrix(table(ancrec.conv$hog, ancrec.conv$ratite==2))
names(ancrec.conv.gene)=c("other", "ratite")
ancrec.conv.gene$ratite.frac = ancrec.conv.gene$ratite / (ancrec.conv.gene$other + ancrec.conv.gene$ratite)
plot(sort(ancrec.conv.gene$ratite.frac), col=ifelse(sort(ancrec.conv.gene$ratite.frac)>0.5, "red", "black"), ylab="Fraction of Convergent Mutations that are Ratite/Ratite")
