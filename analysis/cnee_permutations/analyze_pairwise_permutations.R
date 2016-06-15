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
allPairwiseResults<-list()

for (result.type in result.types) {
  permFile=paste0(permDir, result.type, ".pairwise.results.gz", sep="")
  loadCommand=paste0("gzip -cd ", permFile, sep="")
  
  realData<-read.table(paste0(inputDir, result.type, sep=""), header=T, sep="\t")
  permData<-fread(loadCommand, header=F)
  names(permData)=c("RK", "OK", "OR", "MK", "MR", "MO", "CK", "CR", "CM", "CO")
  realRes<-data.frame(pair=names(permData), count=numeric(10), pval=numeric(10), es=numeric(10))

  realRes[1,2]=sum((realData[,7] + realData[,4]) == 2, na.rm=T)
  realRes[2,2]=sum((realData[,6] + realData[,4]) == 2, na.rm=T)
  realRes[3,2]=sum((realData[,6] + realData[,7]) == 2, na.rm=T)
  realRes[4,2]=sum((realData[,5] + realData[,4]) == 2, na.rm=T)
  realRes[5,2]=sum((realData[,5] + realData[,6]) == 2, na.rm=T)
  realRes[6,2]=sum((realData[,5] + realData[,7]) == 2, na.rm=T)
  realRes[7,2]=sum((realData[,3] + realData[,4]) == 2, na.rm=T)
  realRes[8,2]=sum((realData[,3] + realData[,7]) == 2, na.rm=T)
  realRes[9,2]=sum((realData[,3] + realData[,5]) == 2, na.rm=T)
  realRes[10,2]=sum((realData[,3] + realData[,6]) == 2, na.rm=T)
    
  
  realRes$pval[1] = (1+sum(permData[,RK] >= realRes$count[1]))/(length(permData[,RK])+1)
  realRes$pval[2] = (1+sum(permData[,OK] >= realRes$count[2]))/(length(permData[,OK])+1)
  realRes$pval[3] = (1+sum(permData[,OR] >= realRes$count[3]))/(length(permData[,OR])+1)
  realRes$pval[4] = (1+sum(permData[,MK] >= realRes$count[4]))/(length(permData[,MK])+1)
  realRes$pval[5] = (1+sum(permData[,MR] >= realRes$count[5]))/(length(permData[,MR])+1)
  realRes$pval[6] = (1+sum(permData[,MO] >= realRes$count[6]))/(length(permData[,MO])+1)
  realRes$pval[7] = (1+sum(permData[,CK] >= realRes$count[7]))/(length(permData[,CK])+1)
  realRes$pval[8] = (1+sum(permData[,CR] >= realRes$count[8]))/(length(permData[,CR])+1)
  realRes$pval[9] = (1+sum(permData[,CM] >= realRes$count[9]))/(length(permData[,CM])+1)
  realRes$pval[10] = (1+sum(permData[,CO] >= realRes$count[10]))/(length(permData[,CO])+1)
  
  realRes$es[1]=realRes$count[1]/median(permData[,RK])
  realRes$es[2]=realRes$count[2]/median(permData[,OK])
  realRes$es[3]=realRes$count[3]/median(permData[,OR])
  realRes$es[4]=realRes$count[4]/median(permData[,MK])
  realRes$es[5]=realRes$count[5]/median(permData[,MR])
  realRes$es[6]=realRes$count[6]/median(permData[,MO])
  realRes$es[7]=realRes$count[7]/median(permData[,CK])
  realRes$es[8]=realRes$count[8]/median(permData[,CR])
  realRes$es[9]=realRes$count[9]/median(permData[,CM])
  realRes$es[10]=realRes$count[10]/median(permData[,CO])
  allPairwiseResults[[result.type]]=realRes
  
}
