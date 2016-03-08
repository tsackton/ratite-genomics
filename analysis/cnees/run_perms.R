#run permutations

#get options
args = commandArgs(trailingOnly=TRUE)
#args[1] = is the permutation number

#set number of reps
nreps<-10000
input<-read.table("cnee.for.perm", sep="\t", stringsAsFactors=F, header=T)

clade.res=matrix(nrow=nreps, ncol=10)
colnames(clade.res)=c("ct0.n", "ct1.n", "ct2.n", "ct3.n", "ct4.n", "ct0.b", "ct1.b", "ct2.b", "ct3.b", "ct4.b")
ens.res.br=matrix(nrow=nreps, ncol=length(unique(input$best_ens)))
colnames(ens.res.br)=sort(unique(input$best_ens))
ncbi.res.br=matrix(nrow=nreps, ncol=length(unique(input$best_ncbi)))
colnames(ncbi.res.br)=sort(unique(input$best_ncbi))
ncbi.res.na = ncbi.res.br
ens.res.na = ens.res.br

for (iter in 1:nreps) {
  #sample
  cnee.perm<-as.data.frame(apply(input, 2, sample), stringsAsFactors=F)
  #get clade counts
  clade.ct.br=apply(cnee.perm[,c("rhea.br", "casuar.br", "kiwi.br", "strcam.br")], 1, function(x) sum(x > 0))
  clade.ct.na=apply(cnee.perm[,c("rhea.na", "casuar.na", "kiwi.na", "strcam.na")], 1, function(x) sum(x > 0))
  clades<-data.frame(x0=sum(clade.ct.na==0), x1=sum(clade.ct.na==1), x2=sum(clade.ct.na==2), x3=sum(clade.ct.na==3), x4=sum(clade.ct.na==4), b0=sum(clade.ct.br==0), b1=sum(clade.ct.br==1), b2=sum(clade.ct.br==2), b3=sum(clade.ct.br==3), b4=sum(clade.ct.br==4))
  clades=as.matrix(clades)
  #get gene counts
  genes.ens.b<-table(cnee.perm$accel.broad == 1, cnee.perm$best_ens)[2,]
  genes.ens.n<-table(cnee.perm$accel.narrow == 1, cnee.perm$best_ens)[2,]
  genes.ncbi.b<-table(cnee.perm$accel.broad == 1, cnee.perm$best_ncbi)[2,]
  genes.ncbi.n<-table(cnee.perm$accel.narrow == 1, cnee.perm$best_ncbi)[2,]
  #final results
  clade.res[iter,]=clades
  ncbi.res.na[iter,]=genes.ncbi.n
  ncbi.res.br[iter,]=genes.ncbi.b
  ens.res.na[iter,]=genes.ens.n
  ens.res.br[iter,]=genes.ens.b
}

write.table(clade.res, file=paste0("output/clade.res.perm",args[1]), row.names=F, sep="\t", quote=F)
write.table(ncbi.res.na, file=paste0("output/ncbi.res.narrow.perm",args[1]), row.names=F, sep="\t", quote=F)
write.table(ncbi.res.br, file=paste0("output/ncbi.res.broad.perm",args[1]), row.names=F, sep="\t", quote=F)
write.table(ens.res.na, file=paste0("output/ens.res.narrow.perm",args[1]), row.names=F, sep="\t", quote=F)
write.table(ens.res.br, file=paste0("output/ens.res.broad.perm",args[1]), row.names=F, sep="\t", quote=F)


