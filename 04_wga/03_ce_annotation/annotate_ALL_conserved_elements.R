#This is code to load the conserved elements into R and process the associated bed files
#Takes a path and a bed file and then processes associated files to produce an annotation file
#We focus on the Mitchell et al tree (ce2)
#REQUIRES A NUMBER OF LARGE FILES NOT ON GITHUB (WILL BE AVAILABLE ON DRYAD)

bedfile<-"final_ces_noratite.tree2.bed"
filepath<-"/Volumes/LaCie/Projects/Current/ratites/final/ce_final/"

ce2<-read.table(file=paste0(filepath, bedfile, ".gz", collapse=""), header=F)
lowe<-read.table(file=paste0(filepath, "lowe_cnees.bed.gz", collapse=""), header=F)

#add lengths to each
ce2$length<-ce2$V3-ce2$V2
lowe$length<-lowe$V3-lowe$V2

#plot length distributions
plot(density(ce2$length[ce2$length<1000]), xlim=c(0,300), col="red", lwd=2, main="Conserved element length distribution", xlab="Length", bty="n")
lines(density(lowe$length[lowe$length<1000]), lwd=2, lty="dashed")
legend("topright", legend=c("Sackton et al", "Lowe et al"), lty=c("solid", "dashed"), col=c("red", "black"), lwd=2, bty="n")
abline(v=median(ce2$length), col="red", lty="dotted", lwd=2)
abline(v=median(lowe$length), col="black", lty="dotted", lwd=2)

#length distibutions of elements > 50 bp
plot(density(ce2$length[ce2$length<2000 & ce2$length > 50]), xlim=c(0,500), col="red", lwd=2, main="Conserved element length distribution (50bp cutoff)", xlab="Length", bty="n")
lines(density(lowe$length[lowe$length<2000 & lowe$length > 50]), lwd=2, lty="dashed")
legend("topright", legend=c("Sackton et al", "Lowe et al"), lty=c("solid", "dashed"), col=c("red", "black"), lwd=2, bty="n")
abline(v=median(ce2$length[ce2$length<2000 & ce2$length > 50]), col="red", lty="dotted", lwd=2)
abline(v=median(lowe$length[lowe$length<2000 & lowe$length > 50]), col="black", lty="dotted", lwd=2)

cds.ct<-read.table(file=paste0(filepath, bedfile, ".ct.CDS.gz", collapse=""), header=F, sep="\t")
cnee.ct<-read.table(file=paste0(filepath, bedfile, ".ct.lowe.gz", collapse=""), header=F, sep="")
noncds.ct<-read.table(file=paste0(filepath, bedfile, ".ct.exon.gz", collapse=""), header=F, sep="\t")

#estimate fraction overlapping different categories
prop.table(table(cds.ct$V5 > 0))
prop.table(table(cnee.ct$V5 > 0))
prop.table(table(noncds.ct$V5 > 0))

#load counts and produce annotation
cds.annot<-read.table(file=paste0(filepath, bedfile, ".annot.CDS.gz", collapse=""), header=F, sep="\t")
exon.annot<-read.table(file=paste0(filepath, bedfile, ".annot.exon.gz", collapse=""), header=F, sep="\t")
gene.annot<-read.table(file=paste0(filepath, bedfile, ".annot.gene.gz", collapse=""), header=F, sep="\t")

intersect<-read.table(file=paste0(filepath, "final_ces.intersection.bed.gz", collapse=""), header=F, sep="\t")
intersect$in.intersection=T

#clean up tables
names(gene.annot)[c(4,5)]=c("id", "in.gene")
names(cds.annot)[c(4,5)]=c("id", "in.cds")
names(exon.annot)[c(4,5)]=c("id", "in.exon")

#merge
ce.annot<-merge(cds.annot[,c(4,5)], exon.annot[,c(4,5)], by="id")
ce.annot<-merge(ce.annot, gene.annot[,c(4,5)], by="id")
ce.annot<-merge(ce.annot, ce2[,c("V4", "length")], by.y="V4", by.x="id")
ce.annot<-merge(ce.annot, intersect[,c("V4", "in.intersection")], all.x=T, all.y=F, by.x="id", by.y="V4")
ce.annot$in.intersection[is.na(ce.annot$in.intersection)]=F
ce.annot$class="intergenic"
ce.annot$class[ce.annot$in.gene>0]="genic_non_exonic"
ce.annot$class[ce.annot$in.exon>0]="exonic_non_CDS"
ce.annot$class[ce.annot$in.cds>0]="CDS"

#get closest genes -- will need to collapse duplicate lines
#each CE has four closest genes (ens + ncbi, whole gene vs TSS)
closest.gene.ncbi <- read.table(file=paste0(filepath, bedfile, ".annot.closest_genes_ncbi.gz", collapse=""), header=F, sep="\t", stringsAsFactors=F)
closest.gene.ens <- read.table(file=paste0(filepath, bedfile, ".annot.closest_genes_ens.gz", collapse=""), header=F, sep="\t", stringsAsFactors=F)
closest.tss.ncbi <- read.table(file=paste0(filepath, bedfile, ".annot.closest_TSS_ncbi.gz", collapse=""), header=F, sep="\t", stringsAsFactors=F)
closest.tss.ens <- read.table(file=paste0(filepath, bedfile, ".annot.closest_TSS_ens.gz", collapse=""), header=F, sep="\t", stringsAsFactors=F)

#deal with duplicates
#basic rule: prefer upstream (negative) to downstream
#if there are two 0s then there is no rational basis for deciding which to assign the gene to, so keep it as a duplicate
closest.all.raw<-list(tss_ens=closest.tss.ens, tss_ncbi=closest.tss.ncbi, gene_ens=closest.gene.ens, gene_ncbi=closest.gene.ncbi)
closest.all.clean<-list()

for (i in 1:4) {
  closest<-closest.all.raw[[i]]  
  closest<-closest[order(closest$V4, closest$V11, decreasing=F),]
  closest.remove_dups = subset(closest, (V11==0) | (!duplicated(V4)))
  closest.remove_dups = closest.remove_dups[,c("V4", "V8")]
  names(closest.remove_dups) = c("id", "gene")
  closest.clean = aggregate(gene ~ id, data=closest.remove_dups, paste0, collapse="|")
  names(closest.clean)[2]=names(closest.all.raw)[i]
  closest.clean=closest.clean[order(closest.clean$id),]
  closest.all.clean[[i]]=closest.clean
}
names(closest.all.clean) = names(closest.all.raw)

#merge
closest.gene.final=do.call("cbind", closest.all.clean)

#final cleanup
closest.gene.final=closest.gene.final[,c(1,2,4,6,8)]
names(closest.gene.final)=c("id", "tss_ens", "tss_ncbi", "gene_ens", "gene_ncbi")

#add closest gene and dist to ce.annot
ce.annot<-merge(ce.annot, closest.gene.final, by="id", all.x=T)

#ce on unlocalized scaffolds with no annotations have no adjacent genes
ce.annot$tss_ens[ce.annot$tss_ens=="."]="none"
ce.annot$tss_ncbi[ce.annot$tss_ncbi=="."]="none"
ce.annot$gene_ens[ce.annot$tgene_ens=="."]="none"
ce.annot$gene_ncbi[ce.annot$gene_ncbi=="."]="none"

#get lowe et al cnee overlap, if any
lowe.overlap <- read.table(file=paste0(filepath, bedfile, ".annot.lowe.gz", collapse=""), stringsAsFactors=F)

#clean up lowe et al overlap
lowe.overlap=lowe.overlap[,c("V4", "V8")]
names(lowe.overlap)=c("id", "cnee")
lowe.overlap$cnee[lowe.overlap$cnee == "."] = "none"
lowe.overlap.clean<-aggregate(cnee ~ id, data=lowe.overlap, paste0, collapse="|")

ce.annot<-merge(ce.annot, lowe.overlap.clean, by="id", all.x=T)

#add best ncbi and best ens --> this uses tss for intergenic and gene otherwise
ce.annot$best_ncbi = ce.annot$gene_ncbi
ce.annot$best_ncbi[ce.annot$class=="intergenic"]=ce.annot$tss_ncbi[ce.annot$class=="intergenic"]
ce.annot$best_ens = ce.annot$gene_ens
ce.annot$best_ens[ce.annot$class=="intergenic"]=ce.annot$tss_ens[ce.annot$class=="intergenic"]

#write table
write.table(ce.annot, file=paste0(filepath, bedfile, ".final_annotation.tsv", collapse=""), row.names=F, sep="\t", quote=F)
