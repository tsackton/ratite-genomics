#run permutations

#get options
args = commandArgs(trailingOnly=TRUE)
#args[1] is the permutation number
#args[2] is the input file

#order should be ratite 2 cas 3 kiwi 4 moa 5 ostr 6 rhea 7

#set number of reps
nreps<-10000
input<-read.table(paste0("inputs/",args[2], sep=""), sep="\t", stringsAsFactors=F, header=T, colClasses=c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))

#convergence across ratites
clade.res=matrix(nrow=nreps, ncol=6)
colnames(clade.res)=c("ct0", "ct1", "ct2", "ct3", "ct4", "ct5")

#gene-based permutations
ens.res=matrix(nrow=nreps, ncol=length(unique(input$best_ens)))
colnames(ens.res)=sort(unique(input$best_ens))

#pairwise
pairwise.res=matrix(nrow=nreps, ncol=10)
colnames(pairwise.res)=c("RK", "OK", "OR", "MK", "MR", "MO", "CK", "CR", "CM", "CO")

cnee.perm<-input

for (iter in 1:nreps) {
  #sample
  cnee.perm[,c(2,3,4,5,6,7)]<-apply(cnee.perm[,c(2,3,4,5,6,7)],2,sample)
  
  #get clade counts
  clade.ct=rowSums(cnee.perm[,3:7], na.rm=T)
  clade.res[iter,]=c(sum(clade.ct==0), sum(clade.ct==1), sum(clade.ct==2), sum(clade.ct==3), sum(clade.ct==4), sum(clade.ct==5))
  
  #cas 3 kiwi 4 moa 5 ostr 6 rhea 7
  #get pairwise counts
  pairwise.res[iter,c("RK")]=sum((cnee.perm[,7] + cnee.perm[,4]) == 2, na.rm=T)
  pairwise.res[iter,c("OK")]=sum((cnee.perm[,6] + cnee.perm[,4]) == 2, na.rm=T)
  pairwise.res[iter,c("OR")]=sum((cnee.perm[,6] + cnee.perm[,7]) == 2, na.rm=T)
  pairwise.res[iter,c("MK")]=sum((cnee.perm[,5] + cnee.perm[,4]) == 2, na.rm=T)
  pairwise.res[iter,c("MR")]=sum((cnee.perm[,5] + cnee.perm[,6]) == 2, na.rm=T)
  pairwise.res[iter,c("MO")]=sum((cnee.perm[,5] + cnee.perm[,7]) == 2, na.rm=T)
  pairwise.res[iter,c("CK")]=sum((cnee.perm[,3] + cnee.perm[,4]) == 2, na.rm=T)
  pairwise.res[iter,c("CR")]=sum((cnee.perm[,3] + cnee.perm[,7]) == 2, na.rm=T)
  pairwise.res[iter,c("CM")]=sum((cnee.perm[,3] + cnee.perm[,5]) == 2, na.rm=T)
  pairwise.res[iter,c("CO")]=sum((cnee.perm[,3] + cnee.perm[,6]) == 2, na.rm=T)
  
  #get gene counts
  ens.res[iter,]=table(cnee.perm[,2] == 1, cnee.perm$best_ens)[2,]
}

write.table(clade.res, file=paste0("output/",args[2],"/clade.res.perm",args[1]), row.names=F, sep="\t", quote=F)
write.table(ens.res, file=paste0("output/",args[2],"/ens.res.perm",args[1]), row.names=F, sep="\t", quote=F)
write.table(pairwise.res, file=paste0("output/",args[2],"/pairwise.res.perm",args[1]), row.names=F, sep="\t", quote=F)


