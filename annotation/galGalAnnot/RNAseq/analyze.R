#read data
timepoints<-c("HH18", "HH20", "HH22")
ens<-list()
ncbi<-list()
for (stage in timepoints) {
  ens[[stage]]<-read.table(paste0("/Volumes/LaCie/Projects/Current/ratities/final/galGal_rnaseq/galGal_ens_", stage, "/abundance.tsv"), header=T)
  ncbi[[stage]]<-read.table(paste0("/Volumes/LaCie/Projects/Current/ratities/final/galGal_rnaseq/galGal_ncbi_", stage, "/abundance.tsv"), header=T) 
}

ens.all<-do.call("cbind", ens)
ncbi.all<-do.call("cbind", ncbi)

#get keys
ens.key<-read.table("/Volumes/LaCie/Projects/Current/ratities/final/galGal_rnaseq/ens_trans_gene.txt", header=T, sep="\t")
names(ens.key)=c("gene_id", "transcript_id")

ncbi.key<-read.table("/Volumes/LaCie/Projects/Current/ratities/final/galGal_rnaseq/ncbi_gene_to_trans.txt", header=F)
names(ncbi.key)=c("gene_id", "transcript_id")
ncbi.key=unique(ncbi.key)

#clean up
ens.clean=subset(ens.all, select=c("HH18.target_id", "HH18.length", "HH18.est_counts", "HH18.tpm", "HH20.est_counts", "HH20.tpm", "HH22.est_counts", "HH22.tpm"))
names(ens.clean)[1:2]=c("transcript_id", "length")
ens.clean$transcript_id=sub("\\.\\d+", "", ens.clean$transcript_id)
ens.clean<-merge(ens.clean, ens.key, by="transcript_id", all.x=T)
ens.clean=ens.clean[,c(9,1,2,3,5,7,4,6,8)]

ncbi.clean=subset(ncbi.all, select=c("HH18.target_id", "HH18.length", "HH18.est_counts", "HH18.tpm", "HH20.est_counts", "HH20.tpm", "HH22.est_counts", "HH22.tpm"))
names(ncbi.clean)[1:2]=c("transcript_id", "length")
ncbi.clean$transcript_id = sub("gi\\|\\d+\\|\\w+\\|", "", ncbi.clean$transcript_id, perl=T)
ncbi.clean$transcript_id = sub("\\|", "", ncbi.clean$transcript_id, perl=T)

ncbi.clean<-merge(ncbi.clean, ncbi.key, by="transcript_id", all.x=T)
ncbi.clean=ncbi.clean[,c(9,1,2,3,5,7,4,6,8)]

#make gene keys by summing across transcripts
library(plyr)
ens.gene<-ddply(ens.clean, .(gene_id), summarize, HH18=sum(HH18.tpm), HH20=sum(HH20.tpm), HH22=sum(HH22.tpm))
ncbi.gene<-ddply(ncbi.clean, .(gene_id), summarize, HH18=sum(HH18.tpm), HH20=sum(HH20.tpm), HH22=sum(HH22.tpm))

#read in ce data
conv.cnee<-read.table("convergent_cnees.txt", header=F, sep="\t", fill=T)
names(conv.cnee)=c("gene", "count", "length", "symbol")
conv.cnee<-merge(conv.cnee, ens.gene, by.x="gene", by.y="gene_id", all.x=T, all.y=F)

conv.gene<-read.table("sig_conv_genes.txt", header=F, sep="\t", fill=T)
conv.gene<-merge(conv.gene, ens.gene, by.x="V1", by.y="gene_id", all.x=T, all.y=F)

sig.gene<-read.table("sig_overrep_genes.txt", header=F, sep="\t", fill=T)
sig.gene<-merge(sig.gene, ens.gene, by.x="V1", by.y="gene_id", all.x=T, all.y=F)

#output with expression
sig.gene$type="enriched for accelerated CNEEs at gene level"
conv.gene$type="significant gene level convergence"
conv.cnee$type="significant cnee-level convergence"
names(sig.gene)[1:2]=c("gene", "symbol")
names(conv.gene)[1:2]=c("gene", "symbol")

conv.cnee.sub=unique(conv.cnee[,c(1,4,5,6,7,8)])

all.candidates<-rbind(sig.gene, conv.gene, conv.cnee.sub)

write.table(all.candidates, file="cnee_candidates_with_expression.txt", sep="\t", quote=F, row.names=F)
write.table(ens.gene, file="ens_expression.txt", sep="\t", quote=F, row.names=F)
write.table(ncbi.gene, file="ncbi_expression.txt", sep="\t", quote=F, row.names=F)
