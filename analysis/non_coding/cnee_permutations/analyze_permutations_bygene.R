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

inputDir="./inputs/"
permDir="./"
setwd("./")

permHeader = scan("ens.res.header", what="character")

for (result.type in result.types) {
  permFile=paste0(permDir, result.type, ".gene.results.gz", sep="")
  loadCommand=paste0("gzip -cd ", permFile, sep="")
  
  realData<-read.table(paste0(inputDir, result.type, sep=""), header=T, sep="\t")
  permData<-fread(loadCommand, header=F)
  names(permData) = permHeader

  realData$ct=apply(realData[,2:7], 1, sum, na.rm=T)
  realRes<-as.data.frame(table(realData$ct > 0, realData$best_ens)[2,])
  names(realRes)="ct"
  for (gene in rownames(realRes)) {
    realRes[gene,"pval"] = (1+sum(permData[,gene,with=FALSE] >= realRes[gene,"ct"]))/(nrow(permData[,gene,with=FALSE])+1)
  }  
  write.table(realRes, file=paste0(result.type, ".gene_pvals", sep=""), sep="\t", quote=F)
}
