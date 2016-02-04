##This is code used to test the sensitivity of different phast parameters##
##Not used to analyze final dataset##

#make target coverage / expected length data frame for each iteration

meta<-data.frame(iteration=seq(1,9,1), targetcov=c(0.3,0.4,0.4,0.4,0.2,0.2,0.2,0.3,0.3), explen=c(45,20,45,70,20,45,70,70,20))

setwd("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast")

beds<-list()
gene.cov<-list()
cnee.cov<-list()
for (ite in meta$iteration) {
  bedfile = paste("phastCons_output/", "cons_ite", ite, ".bed", sep="")
  genefile = paste("phastCons_output/", "ite", ite, ".gene.count", sep="")
  cneefile = paste("phastCons_output/", "ite", ite, ".lowe_cnee.count", sep="")
  print(bedfile)
 beds[[ite]]<-read.delim(bedfile, stringsAsFactors=F, header=F)
 gene.cov[[ite]]<-read.delim(genefile, stringsAsFactors=F, header=F)
 cnee.cov[[ite]]<-read.delim(cneefile, stringsAsFactors=F, header=F)
}

lowe.cnee<-read.delim("phastCons_output/lowe_cnees.bed", stringsAsFactors=F, header=F)
lowe.cnee$length <- lowe.cnee$V3 - lowe.cnee$V2
chicken.exons<-read.delim("phastCons_output/chicken_genes.bed", stringsAsFactors=F, header=F)
chicken.exons$length <- chicken.exons$V3 - chicken.exons$V2

for (ite in meta$iteration) {
  beds[[ite]]$length = beds[[ite]]$V3 - beds[[ite]]$V2
}

pdf(file="length_dists.pdf")
max.len<-500
lengths<-lowe.cnee$length
lengths[lengths > max.len]=max.len
hist(lengths, breaks=50, main="Lowe CNEEs", xlab="Length")
num.gt.100 = sum(lengths>100)
legend("topright", legend=c(paste0("Num > 100: ", num.gt.100)), lwd=0, col="white", bty="n")
lengths<-chicken.exons$length
lengths[lengths > max.len]=max.len
hist(lengths, breaks=50, main="Chicken Exons", xlab="Length")
num.gt.100 = sum(lengths>100)
legend("topright", legend=c(paste0("Num > 100: ", num.gt.100)), lwd=0, col="white", bty="n")

for (i in meta[order(meta$targetcov, meta$explen), "iteration"]) {
    tc<-meta$targetcov[i]
    ite<-meta$iteration[i]
    el<-meta$explen[i]
    lengths<-beds[[i]]$length
    lengths[lengths > max.len]=max.len
    hist(lengths, breaks=50, main=paste("Target Coverage: ", tc, " Expected Length: ", el, sep=""), xlab="Length")
    num.gt.100 = sum(lengths>100)
    legend("topright", legend=c(paste0("Num > 100: ", num.gt.100)), lwd=0, col="white", bty="n")
}
dev.off()

pdf(file="cons_per_exon.pdf")
max.genecount<-25
for (i in meta[order(meta$targetcov, meta$explen), "iteration"]) {
  tc<-meta$targetcov[i]
  ite<-meta$iteration[i]
  el<-meta$explen[i]
  genecounts<-gene.cov[[ite]]$V7
  genecounts[genecounts > max.genecount]=max.genecount
  hist(genecounts, breaks=seq(0,max.genecount+1,1), include.lowest=F, right=F, xlab="Number of conserved elements per exon", main=paste("Target Coverage: ", tc, " Expected Length: ", el, sep=""))
  num.zero = round(sum(genecounts==0)/length(genecounts),4)
  num.gt.2 = round(sum(genecounts>=2)/length(genecounts),4)
  num.one = round(1-num.zero-num.gt.2,4)
  legend("topright", legend=c(paste0("Frac == 0: ", num.zero), paste0("Frac == 1: ", num.one), paste0("Frac >= 2: ", num.gt.2)), lwd=0, col="white", bty="n")
}
dev.off()


pdf(file="cons_per_cnee.pdf")
max.genecount<-10
for (i in meta[order(meta$targetcov, meta$explen), "iteration"]) {
  tc<-meta$targetcov[i]
  ite<-meta$iteration[i]
  el<-meta$explen[i]
  genecounts<-cnee.cov[[ite]]$V5
  genecounts[genecounts > max.genecount]=max.genecount
  hist(genecounts, breaks=seq(0,max.genecount+1,1), include.lowest=F, right=F, xlab="Number of conserved elements per exon", main=paste("Target Coverage: ", tc, " Expected Length: ", el, sep=""))
  num.zero = round(sum(genecounts==0)/length(genecounts),4)
  num.gt.2 = round(sum(genecounts>=2)/length(genecounts),4)
  num.one = round(1-num.zero-num.gt.2,4)
  legend("topright", legend=c(paste0("Frac == 0: ", num.zero), paste0("Frac == 1: ", num.one), paste0("Frac >= 2: ", num.gt.2)), lwd=0, col="white", bty="n")
}
dev.off()




