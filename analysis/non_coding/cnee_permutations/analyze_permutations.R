#read inputs
library(data.table)

result.types<-character(0)

for (align in c("with_moa", "wga")) {
  for (model in c("neut_ver1", "neut_ver2", "neut_ver3")) {
    for (type in c("1", "2")) {
      result.types=c(result.types, paste0(align, ".", model, ".for_perm.", type))
    }
  }
}

inputDir="~/Projects/birds/ratite_compgen/ratite-genomics/analysis/cnee_permutations/inputs/"
permDir="/Volumes/LaCie/Projects/Current/ratites/final/cnee_perms/"
setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/cnee_permutations/")

for (result.type in result.types) {
  permFile=paste0(permDir, result.type, ".clade.results.gz", sep="")
  loadCommand=paste0("gzip -cd ", permFile, sep="")
  
  realData<-read.table(paste0(inputDir, result.type, sep=""), header=T, sep="\t")
  permData<-fread(loadCommand, header=F)
  realData$ct=apply(realData[,3:7], 1, sum, na.rm=T)
  realRes<-data.frame(ct=c(1,2,3,4,5), num=numeric(5))
  
  realRes$num[1]=sum(realData$ct >= 1)
  realRes$num[2]=sum(realData$ct >= 2)
  realRes$num[3]=sum(realData$ct >= 3)
  realRes$num[4]=sum(realData$ct >= 4)
  realRes$num[5]=sum(realData$ct >= 5)
  
  #get pvalues
  permData$ct1 = permData$V2+permData$V3+permData$V4+permData$V5+permData$V6
  permData$ct2 = permData$V3+permData$V4+permData$V5+permData$V6
  permData$ct3 = permData$V4+permData$V5+permData$V6
  permData$ct4 = permData$V5+permData$V6
  permData$ct5 = permData$V6
  
  realRes$pval[1] = (1+sum(permData[,ct1] >= realRes$num[1]))/(length(permData[,ct1])+1)
  realRes$pval[2] = (1+sum(permData[,ct2] >= realRes$num[2]))/(length(permData[,ct2])+1)
  realRes$pval[3] = (1+sum(permData[,ct3] >= realRes$num[3]))/(length(permData[,ct3])+1)
  realRes$pval[4] = (1+sum(permData[,ct4] >= realRes$num[4]))/(length(permData[,ct4])+1)
  realRes$pval[5] = (1+sum(permData[,ct5] >= realRes$num[5]))/(length(permData[,ct5])+1)
  
  #plot histograms
  pdf(file=paste0(result.type, ".clade.pdf", sep=""))
  for (index in c(1, 2, 3, 4, 5)) {
    xmin=min(unlist(c(permData[,paste0("ct",index,sep=""),with=F], realRes$num[index])))*0.97
    xmax=max(unlist(c(permData[,paste0("ct",index,sep=""),with=F], realRes$num[index])))*1.03
    breaks=min(100,nrow((unique(permData[,paste0("ct",index,sep=""),with=F]))))
    
    hist(unlist(permData[,paste0("ct",index,sep=""),with=F]), breaks=breaks, xlim=c(xmin,xmax), xlab="Values", main=paste0("Clade Count ", realRes$ct[index], " [P(obs>null) = ", realRes$pval[index], "]", sep=""), col="blue", border="blue")
    abline(v=realRes$num[index], lwd=3, col="red")
  }  
  dev.off()  
}
