#CNEE analysis, final results for manuscript

#load data -- this is rather slow so should refactor to save actual datasets and just load what I need, I guess
source("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/cnees/analyze_cnees.R")
source("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/cnee_permutations/analyze_permutations.R")
source("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/cnees/generate_candidates.R")

#linear model testing of probability of loss in each ratite lineage
#use accel.clade, but subset to be just with_moa.neut_ver3
accel.moa3<-subset(accel.clade, result.type=="with_moa.neut_ver3")

#logit regression test
ostrich.logit<-glm(ostrich.acc.2 ~ nonAvesPres + AvesPres + cas.acc.2 + kiwi.acc.2 + rhea.acc.2 + class + length + dupcount, data=accel.moa3, family="binomial")
kiwi.logit<-glm(kiwi.acc.2 ~ nonAvesPres +  AvesPres + cas.acc.2 + ostrich.acc.2 + rhea.acc.2 + class + length + dupcount, data=accel.moa3, family="binomial")
rhea.logit<-glm(rhea.acc.2 ~ nonAvesPres +  AvesPres + cas.acc.2 + kiwi.acc.2 + ostrich.acc.2 + class + length + dupcount, data=accel.moa3, family="binomial")
cas.logit<-glm(cas.acc.2 ~ nonAvesPres +  AvesPres + ostrich.acc.2 + kiwi.acc.2 + rhea.acc.2 + class + length + dupcount, data=accel.moa3, family="binomial")
summary(ostrich.logit)
summary(kiwi.logit)
summary(rhea.logit)
summary(cas.logit)

ratite.logit.res<-data.frame(species=rep(c("ostrich", "rhea", "kiwi", "emu/cassowary"), each=3), est=numeric(12), stderr=numeric(12))
ratite.logit.res[1:3,2:3]=coef(summary(ostrich.logit))[4:6,c(1,2)]
ratite.logit.res[4:6,2:3]=coef(summary(rhea.logit))[4:6,c(1,2)]
ratite.logit.res[7:9,2:3]=coef(summary(kiwi.logit))[4:6,c(1,2)]
ratite.logit.res[10:12,2:3]=coef(summary(cas.logit))[4:6,c(1,2)]

ratite.logit.res$ci.u = ratite.logit.res$est + (1.96*ratite.logit.res$stderr)
ratite.logit.res$ci.l = ratite.logit.res$est - (1.96*ratite.logit.res$stderr)
library(plotrix)
plotCI(x=seq(1,12), y=exp(ratite.logit.res$est), ui=exp(ratite.logit.res$ci.u), li=exp(ratite.logit.res$ci.l), xaxt="n", xlab="", bty="n", las=2, ylab="Estimate (logit model)", ylim=c(0,10))
abline(h=1, col="red", lty="dashed")
axis(1, at=c(2, 5, 8, 11), labels=unique(ratite.logit.res$species))
abline(v=3.5, lty="dashed", col="gray")
abline(v=6.5, lty="dashed", col="gray")
abline(v=9.5, lty="dashed", col="gray")


#gene pvals
genepval<-read.table("cnee_permutations/with_moa.neut_ver2.for_perm.2.gene_pvals")
genetot<-as.data.frame(table(cnees$best_ens))
head(genepval)
genepval<-merge(genepval, genetot, by.x="row.names", by.y="Var1")
genepval$qval<-p.adjust(genepval$pval, method="fdr")
hist(genepval$pval[genepval$pval<0.99], breaks=50, col="red")

plot(genepval$Freq ~ genepval$ct, pch=16, col=ifelse(genepval$qval < 0.1, "red", "black"), log="xy", cex=0.75, las=1, xlab="# ratite-accelerated CNEEs", ylab="Total # CNEEs", bty="n")
