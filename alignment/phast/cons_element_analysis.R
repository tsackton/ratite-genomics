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

#update with emu info, and merge ratite acceleration
ce.annot<-read.table("~/Dropbox/Public/ce_annotation.tsv.gz", header=T)
ratite.accel<-read.table("~/Dropbox/Public/ratite_accel.tsv.gz", header=F)
names(ratite.accel)=c("id", "ratite.pval")
emu.accel<-read.table("final_beds/emu_accel_sorted.out", header=F)
emu.accel=emu.accel[,c("V4", "V9")]
names(emu.accel)=c("id", "emu.pval")

ce.merge<-merge(ce.annot, ratite.accel, by.x="id", by.y="id")
ce.merge<-merge(ce.merge, emu.accel, by="id")

ce.merge$ratite.sign = sign(ce.merge$ratite.pval)
ce.merge$emu.sign = sign(ce.merge$emu.pval)
ce.merge$ratite.q = p.adjust(abs(ce.merge$ratite.pval), method="fdr")
ce.merge$emu.q = p.adjust(abs(ce.merge$emu.pval), method="fdr")
ce.merge$emu.sig = as.numeric(ce.merge$emu.q < 0.05) * ce.merge$emu.sign
ce.merge$ratite.sig = as.numeric(ce.merge$ratite.q < 0.05) * ce.merge$ratite.sign

head(ce.merge)

#now keep track at the individual element level
ce.merge$element.key = "none"
ce.merge$element.key[ce.merge$emu.sig < 0 & ce.merge$ratite.sig < 0]="e_and_r"
ce.merge$element.key[ce.merge$emu.sig < 0 & ce.merge$ratite.sig >= 0]="e_only"
ce.merge$element.key[ce.merge$emu.sig >= 0 & ce.merge$ratite.sig < 0]="r_only"

cnee <- subset(ce.merge, (class == "genic_non_exonic" | class == "intergenic") & length > 50)

sig.table<-table(cnee$gene, cnee$element.key)
sig.df<-as.data.frame.matrix(sig.table)
sig.df$total=sig.df$e_and_r+sig.df$e_only+sig.df$r_only+sig.df$none

ncbi.to.ens<-read.table("~/Downloads/CGNC_Gallus_gallus_20151214.txt", header=F, sep="\t", quote="", comment.char="")
names(ncbi.to.ens)=c("cgnc", "ncbi", "ncbi.ver", "ens", "sym", "name", "specie")
ncbi.to.ens=subset(ncbi.to.ens, select=c("cgnc", "ncbi", "ens", "sym"))

sig.df$cgnc=sub(".*CGNC:(\\d+).*", "\\1", rownames(sig.df), perl=T)
sig.df$ncbi=sub(".*GeneID:(\\d+).*", "\\1", rownames(sig.df), perl=T)

sig.df$cgnc[grepl("GeneID", sig.df$cgnc)]=NA
sig.df$ncbi[grepl("CGNC", sig.df$ncbi)]=NA

sig.df=sig.df[-1,]

#add ens and symbol
sig.df=merge(sig.df, ncbi.to.ens, by="ncbi", all.x=T, all.y=F)
sig.df$ens[sig.df$ens==""]=NA
sig.df$cgnc=sig.df$cgnc.y
sig.df=sig.df[,c(1:6,9,10,11)]
head(sig.df)

#duplicate rows happen because of the parsing of the 'ties'
#just sum them
#do this by aggregating by ncbi ID then merging
annotation<-sig.df[,c("ncbi", "ens", "sym", "cgnc")]
annotation<-unique(annotation)
data<-sig.df[,c("ncbi", "e_and_r", "e_only", "r_only", "none", "total")]
library(plyr)
merged.data<-ddply(data, "ncbi", numcolwise(sum))

final.data<-merge(merged.data, annotation, by.x="ncbi", by.y="ncbi", all.x=T, all.y=F)

#a few duplicates remain due to som disagreements about sym/cgnc id
#remove by hand even though it is ugly
final.data[final.data$ncbi %in% final.data$ncbi[duplicated(final.data$ncbi)],]
final.data=final.data[-16187,]
final.data=final.data[-14171,]
final.data=final.data[-4825,]
final.data=final.data[-853,]
final.data=final.data[-637,]

#check
final.data[final.data$ncbi %in% final.data$ncbi[duplicated(final.data$ncbi)],]

#now add expression
exp<-read.table("~/Downloads/limb_expression.txt", header=T)
head(exp)
final.data.exp<-merge(final.data, exp, by="ens", all.x=T, all.y=F)
final.data.exp$total.accel=final.data.exp$e_and_r+final.data.exp$e_only+final.data.exp$r_only
final.data.exp$avg.exp = apply(final.data.exp[,c("HH18_WT", "HH20_WT", "HH22_WT")], 1, mean, na.rm=T)

table(is.finite(final.data.exp$avg.exp), final.data.exp$total.accel)

write.table(file="ratite_data_ver1.tsv", final.data.exp, row.names=F, quote=F, sep="\t")
