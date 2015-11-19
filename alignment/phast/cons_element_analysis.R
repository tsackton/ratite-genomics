#Analyze conserved elements

ce1<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/most_conserved_final.tree1.bed", header=F)
ce2<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/most_conserved_final.tree2.bed", header=F)
lowe<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/lowe_cnees.bed", header=F)

#add lengths to each
ce1$length<-ce1$V3-ce1$V2
ce2$length<-ce2$V3-ce2$V2
lowe$length<-lowe$V3-lowe$V2

#plot length distributions
plot(density(ce2$length[ce2$length<1000]), xlim=c(0,300), col="red", lwd=2, main="Conserved element length distribution", xlab="Length", bty="n")
lines(density(lowe$length[lowe$length<1000]), lwd=2, lty="dashed")
legend("topright", legend=c("Sackton et al", "Lowe et al"), lty=c("solid", "dashed"), col=c("red", "black"), lwd=2, bty="n")
abline(v=median(ce2$length), col="red", lty="dotted", lwd=2)
abline(v=median(lowe$length), col="black", lty="dotted", lwd=2)

#length distibutions of elements > 50 bp
plot(density(ce2$length[ce2$length<1000 & ce2$length > 75]), xlim=c(0,300), col="red", lwd=2, main="Conserved element length distribution (75bp cutoff)", xlab="Length", bty="n")
lines(density(lowe$length[lowe$length<1000 & lowe$length > 75]), lwd=2, lty="dashed")
legend("topright", legend=c("Sackton et al", "Lowe et al"), lty=c("solid", "dashed"), col=c("red", "black"), lwd=2, bty="n")
abline(v=median(ce2$length[ce2$length<1000 & ce2$length > 75]), col="red", lty="dotted", lwd=2)
abline(v=median(lowe$length[lowe$length<1000 & lowe$length > 75]), col="black", lty="dotted", lwd=2)

cds.ct<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/tree2.CDS.count", header=F, sep="\t")
cnee.ct<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/tree2.lowe_cnee.count", header=F, sep="\t")
noncds.ct<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/tree2.nonCDS.count", header=F, sep="\t")

prop.table(table(cds.ct$V5 > 0))
prop.table(table(cnee.ct$V5 > 0))
prop.table(table(noncds.ct$V5 > 0))

#load counts and produce annotation
cds.annot<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/CDS.annot", header=F, sep="\t")
exon.annot<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/exon.annot", header=F, sep="\t")
gene.annot<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/gene.annot", header=F, sep="\t")

intersect<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/most_conserved_final.intersection.bed", header=F, sep="\t")
intersect$in.intersection=T

#clean up tables
names(gene.annot)[c(4,7)]=c("id", "in.gene")
names(cds.annot)[c(4,7)]=c("id", "in.cds")
names(exon.annot)[c(4,7)]=c("id", "in.exon")

#merge
ce.annot<-merge(cds.annot[,c(4,7)], exon.annot[,c(4,7)], by="id")
ce.annot<-merge(ce.annot, gene.annot[,c(4,7)], by="id")
ce.annot<-merge(ce.annot, ce2[,c("V4", "length")], by.y="V4", by.x="id")
ce.annot<-merge(ce.annot, intersect[,c("V4", "in.intersection")], all.x=T, by.x="id", by.y="V4")
ce.annot$in.intersection[is.na(ce.annot$in.intersection)]=F
ce.annot$class="intergenic"
ce.annot$class[ce.annot$in.gene>0]="genic_non_exonic"
ce.annot$class[ce.annot$in.exon>0]="exonic_non_CDS"
ce.annot$class[ce.annot$in.cds>0]="CDS"
ce.annot<-merge(ce.annot, ce2[,c("V4", "V5")], by.x="id", by.y="V4")

names(ce.annot)[8]="score"

#get closest genes -- will need to collapse duplicate lines
closest.gene <- read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/tree2.closest_genes.out", header=F, sep="\t")
closest.gene<-closest.gene[,c("V4", "V10", "V13")]
names(closest.gene)<-c("id", "gene", "dist")
closest.gene$gene<-sub(",Start\\S+", "", closest.gene$gene, perl=T)  
#collapse duplicates
closest.clean<-aggregate(gene ~ id, data=closest.gene, paste0, collapse="|")
closest.clean<-merge(closest.clean, closest.gene[,c("id", "dist")], all.x=T, by="id")
#some duplicates are equidistant from one gene upstream and one gene downstream,
closest.clean$dist[duplicated(closest.clean$id)]=-1*closest.clean$dist[duplicated(closest.clean$id)]
closest.clean=unique(closest.clean)

#add closest gene and dist to ce.annot
ce.annot<-merge(ce.annot, closest.clean, by="id", all.x=T)

#get lowe et al cnee overlap, if any
lowe.overlap <- read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/lowe.annot", stringsAsFactors=F)

#clean up lowe et al overlap
lowe.overlap=lowe.overlap[,c("V4", "V10")]
names(lowe.overlap)=c("id", "cnee")
lowe.overlap$cnee[lowe.overlap$cnee == "."] = "none"
lowe.overlap.clean<-aggregate(cnee ~ id, data=lowe.overlap, paste0, collapse="|")

ce.annot<-merge(ce.annot, lowe.overlap.clean, by="id", all.x=T)

#write table
write.table(ce.annot, file="ce_annotation.tsv", row.names=F, sep="\t", quote=F)
