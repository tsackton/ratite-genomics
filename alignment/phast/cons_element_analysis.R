#Analyze conserved elements

ce1<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/final_ces.tree1.bed.gz", header=F)
ce2<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/final_ces.tree2.bed.gz", header=F)
ceR<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/final_ces_noratite.tree2.bed.gz", header=F)
lowe<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/lowe_cnees.bed.gz", header=F)

#add lengths to each
ce1$length<-ce1$V3-ce1$V2
ce2$length<-ce2$V3-ce2$V2
ceR$length<-ceR$V3-ceR$V2
lowe$length<-lowe$V3-lowe$V2

#plot length distributions
plot(density(ce2$length[ce2$length<1000]), xlim=c(0,300), col="red", lwd=2, main="Conserved element length distribution", xlab="Length", bty="n")
lines(density(lowe$length[lowe$length<1000]), lwd=2, lty="dashed")
legend("topright", legend=c("Sackton et al", "Lowe et al"), lty=c("solid", "dashed"), col=c("red", "black"), lwd=2, bty="n")
abline(v=median(ce2$length), col="red", lty="dotted", lwd=2)
abline(v=median(lowe$length), col="black", lty="dotted", lwd=2)

#length distibutions of elements > 75 bp
plot(density(ce2$length[ce2$length<2000 & ce2$length > 50]), xlim=c(0,500), col="red", lwd=2, main="Conserved element length distribution (50bp cutoff)", xlab="Length", bty="n")
lines(density(lowe$length[lowe$length<2000 & lowe$length > 50]), lwd=2, lty="dashed")
legend("topright", legend=c("Sackton et al", "Lowe et al"), lty=c("solid", "dashed"), col=c("red", "black"), lwd=2, bty="n")
abline(v=median(ce2$length[ce2$length<2000 & ce2$length > 50]), col="red", lty="dotted", lwd=2)
abline(v=median(lowe$length[lowe$length<2000 & lowe$length > 50]), col="black", lty="dotted", lwd=2)

cds.ct<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/final_ces.tree2.bed.ct.CDS.gz", header=F, sep="\t")
cnee.ct<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/final_ces.tree2.bed.ct.lowe.gz", header=F, sep="")
noncds.ct<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/final_ces.tree2.bed.ct.exon.gz", header=F, sep="\t")

prop.table(table(cds.ct$V5 > 0))
prop.table(table(cnee.ct$V5 > 0))
prop.table(table(noncds.ct$V5 > 0))

#load counts and produce annotation
cds.annot<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/final_ces.tree2.bed.annot.CDS.gz", header=F, sep="\t")
exon.annot<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/final_ces.tree2.bed.annot.exon.gz", header=F, sep="\t")
gene.annot<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/final_ces.tree2.bed.annot.gene.gz", header=F, sep="\t")

intersect<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/final_ces.intersection.bed.gz", header=F, sep="\t")
intersect$in.intersection=T

#clean up tables
names(gene.annot)[c(4,5)]=c("id", "in.gene")
names(cds.annot)[c(4,5)]=c("id", "in.cds")
names(exon.annot)[c(4,5)]=c("id", "in.exon")

#merge
ce.annot<-merge(cds.annot[,c(4,5)], exon.annot[,c(4,5)], by="id")
ce.annot<-merge(ce.annot, gene.annot[,c(4,5)], by="id")
ce.annot<-merge(ce.annot, ce2[,c("V4", "length")], by.y="V4", by.x="id")
ce.annot<-merge(ce.annot, intersect[,c("V4", "in.intersection")], all.x=T, by.x="id", by.y="V4")
ce.annot$in.intersection[is.na(ce.annot$in.intersection)]=F
ce.annot$class="intergenic"
ce.annot$class[ce.annot$in.gene>0]="genic_non_exonic"
ce.annot$class[ce.annot$in.exon>0]="exonic_non_CDS"
ce.annot$class[ce.annot$in.cds>0]="CDS"

#get closest genes -- will need to collapse duplicate lines
closest.gene.ncbi <- read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/final_ces.tree2.bed.annot.closest_genes_ncbi.gz", header=F, sep="\t", stringsAsFactors=F)
closest.gene.ens <- read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/final_ces.tree2.bed.annot.closest_genes_ens.gz", header=F, sep="\t", stringsAsFactors=F)

#run with NCBI
closest.gene<-closest.gene.ncbi[,c("V4", "V8", "V11")]
names(closest.gene)<-c("id", "gene", "dist")
#collapse duplicates
closest.clean<-aggregate(gene ~ id, data=closest.gene, paste0, collapse="|")
closest.clean<-merge(closest.clean, closest.gene[,c("id", "dist")], all.x=T, by="id")
#some duplicates are equidistant from one gene upstream and one gene downstream,
closest.clean$dist[duplicated(closest.clean$id)]=-1*closest.clean$dist[duplicated(closest.clean$id)]
closest.clean.ncbi=unique(closest.clean)

#run with ens
closest.gene<-closest.gene.ens[,c("V4", "V8", "V11")]
names(closest.gene)<-c("id", "gene", "dist")
#collapse duplicates
closest.clean<-aggregate(gene ~ id, data=closest.gene, paste0, collapse="|")
closest.clean<-merge(closest.clean, closest.gene[,c("id", "dist")], all.x=T, by="id")
#some duplicates are equidistant from one gene upstream and one gene downstream,
closest.clean$dist[duplicated(closest.clean$id)]=-1*closest.clean$dist[duplicated(closest.clean$id)]
closest.clean.ens=unique(closest.clean)

#merge
names(closest.clean.ens)=c("id", "ens", "ens.dist")
names(closest.clean.ncbi)=c("id", "ncbi", "ncbi.dist")
closest.gene.final<-merge(closest.clean.ens, closest.clean.ncbi, by="id", all=T)

#add closest gene and dist to ce.annot
ce.annot<-merge(ce.annot, closest.gene.final, by="id", all.x=T)

#get lowe et al cnee overlap, if any
lowe.overlap <- read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/final_ces.tree2.bed.annot.lowe.gz", stringsAsFactors=F)

#clean up lowe et al overlap
lowe.overlap=lowe.overlap[,c("V4", "V8")]
names(lowe.overlap)=c("id", "cnee")
lowe.overlap$cnee[lowe.overlap$cnee == "."] = "none"
lowe.overlap.clean<-aggregate(cnee ~ id, data=lowe.overlap, paste0, collapse="|")

ce.annot<-merge(ce.annot, lowe.overlap.clean, by="id", all.x=T)

#write table
write.table(ce.annot, file="./final_beds/ce_annotation.tsv", row.names=F, sep="\t", quote=F)

#process acceleration results

#load info
accel.classes<-c("Casuar", "Rhea", "Kiwi", "tinamou", "ratite", "aptOwe", "aptHaa", "aptRow", "casCas", "droNov", "rheAme", "rhePen", "strCam")
accel.res<-list()
for (group in accel.classes) {
  accel.res[[group]]<-read.table(paste0("final_accel/all_", group, ".out", sep=""), header=T, sep="\t", comment.char="")
}

for (group in accel.classes) {
  accel.res[[group]]$qval = p.adjust(accel.res[[group]]$pval, method="fdr")
  accel.res[[group]]$group = group
}

accel.all<-do.call("rbind", accel.res)
accel.sub<-subset(accel.all, select=c("name", "qval", "group"))
library(tidyr)
accel.wide<-spread(accel.sub, group, qval)

#make tinamou filter
accel.wide$tin.filt = 0
accel.wide$tin.filt[accel.wide$tinamou <= 0.1] = 1

#get ratite species count
accel.wide$sp.count = apply(accel.wide[,c("aptHaa", "aptOwe", "aptRow", "casCas", "droNov", "strCam", "rheAme", "rhePen")], 1, function(x) sum(x < 0.05))

#merge with ce.annot
ce.merge<-merge(ce.annot, accel.wide, by.x="id", by.y="name")

#subset to remove exonic
cnee<-subset(ce.merge, class == "intergenic" | class == "genic_non_exonic")

cnee$clade.ct = apply(cnee[,c("Rhea", "Casuar", "Kiwi", "strCam")], 1, function(x) sum(x < 0.05))
cnee$total.accel =apply(cnee[,c("aptHaa", "aptOwe", "aptRow", "casCas", "droNov", "strCam", "rheAme", "rhePen", "Rhea", "Casuar", "Kiwi", "ratite")], 1, function(x) sum(x < 0.05))
cnee$ratite.broad = 0
cnee$ratite.narrow = 0
cnee$ratite.broad[cnee$total.accel > 0 & cnee$tin.filt == 0] = 1
cnee$ratite.narrow[cnee$ratite < 0.05 & cnee$tin.filt == 0] = 1

#length subset
cnee.long = subset(cnee, length > 50)

#get elements present in each species of palaeognath
psp <- c("strCam", "droNov", "aptHaa", "aptOwe", "aptRow", "casCas", "rheAme", "rhePen", "tinGut", "cryCin", "eudEle", "notPer")
cnee.pres<-list()
for (group in psp) {
  cnee.pres[[group]]<-read.table(paste0("final_accel/", group, ".bed", sep=""), header=F, sep="\t", comment.char="")
  cnee.pres[[group]]$sp = group
}

cnee.long$tinGutPres = cnee.long$id %in% cnee.pres$tinGut$V4
cnee.long$eudElePres = cnee.long$id %in% cnee.pres$eudEle$V4
cnee.long$notPerPres = cnee.long$id %in% cnee.pres$notPer$V4
cnee.long$cryCinPres = cnee.long$id %in% cnee.pres$cryCin$V4
cnee.long$tinPres = apply(cnee.long[,c("tinGutPres", "eudElePres", "notPerPres", "cryCinPres")], 1, sum)

#top candidates
cnee.long[cnee.long$clade.ct > 2 & cnee.long$ratite.broad == 1 & cnee.long$tinPres > 0,c("id", "ens", "clade.ct", "length")]

cnee.long$ratite.broad.strict = cnee.long$ratite.broad
cnee.long$ratite.broad.strict[cnee.long$tinPres == 0] = 0

#pool across genes -- using ensembl for now
sig.table<-table(cnee.long$ens, cnee.long$ratite.broad.strict)
sig.df<-as.data.frame.matrix(sig.table)
names(sig.df)=c("cons", "accel")
sig.df$total = sig.df$cons + sig.df$accel
sig.df = sig.df[-1,] #remove row with no gene annotation
sig.df$frac = sig.df$accel / sig.df$total
sig.df=sig.df[order(sig.df$frac, decreasing=T),]

#genes with > 10 cnees and > 0.2 accelerated
sig.df[sig.df$total>=10 & sig.df$frac >= 0.2,]

#long list of all genes
write.table(sig.df[sig.df$total>0 & sig.df$accel > 0,], file="all_cnee_genes.tsv", sep="\t", quote=F)

#permutations - version 1, this shuffles the tinamou values as well as the other values
#first make subset that has gene, q-values, and tinPres columns
cnee.input<-subset(cnee.long, select=c("ens", "aptHaa", "aptOwe", "aptRow", "casCas", "Casuar", "droNov", "Kiwi", "ratite", "Rhea", "rheAme", "rhePen", "strCam", "tinamou", "tinPres"))

cnee.input$accel = apply(cnee.input[,c(2:12)], 1, min)
cnee.input=cnee.input[,c("ens", "accel", "Casuar", "Kiwi", "Rhea", "strCam", "tinamou", "tinPres")]

nreps<-100
clade.res=matrix(nrow=nreps, ncol=5)
colnames(clade.res)=c("ct0", "ct1", "ct2", "ct3", "ct4")
gene.res=matrix(nrow=nreps, ncol=length(unique(cnee.input$ens)))
colnames(gene.res)=unique(cnee.input$ens)
genecov.res=gene.res

for (iter in 1:nreps) {
  #shuffle
  cnee.perm<-as.data.frame(apply(cnee.input, 2, sample), stringsAsFactors=F)
  #get clade counts
  clade.ct=apply(cnee.perm[,c("Rhea", "Casuar", "Kiwi", "strCam")], 1, function(x) sum(x < 0.05))
  clade.ct[cnee.perm$tinamou <= 0.1 | cnee.perm$tinPres < 1]=0;
  clades<-data.frame(x0=sum(clade.ct==0), x1=sum(clade.ct==1), x2=sum(clade.ct==2), x3=sum(clade.ct==3), x4=sum(clade.ct==4))
  clades=as.matrix(clades)
  #get gene counts
  genes<-table(cnee.perm$accel < 0.05 & cnee.perm$tinPres > 0 & cnee.perm$tinamou > 0.1, cnee.perm$ens)[2,]
  #gene level convergence (may be slow)
  cnee.perm[cnee.perm$tinamou <=0.1 | cnee.perm$tinPres == 0,c(3,4,5,6)]=1
  gene.conv=ddply(cnee.perm, .(ens), summarize, r.min=min(Rhea), c.min=min(Casuar), o.min=min(strCam), k.min=min(Kiwi))
  gene.conv.ct=apply(gene.conv[,c("r.min", "c.min", "k.min", "o.min")], 1, function(x) sum(x < 0.05))
  names(gene.conv.ct)=gene.conv$ens  
  #final results
  clade.res[iter,]=clades
  gene.res[iter,]=genes
  genecov.res[iter,]=gene.conv.ct
}

#real results
cnee.perm<-cnee.input
#get clade counts
clade.ct=apply(cnee.perm[,c("Rhea", "Casuar", "Kiwi", "strCam")], 1, function(x) sum(x < 0.05))
clade.ct[cnee.perm$tinamou <= 0.1 | cnee.perm$tinPres < 1]=0;
clades<-data.frame(x0=sum(clade.ct==0), x1=sum(clade.ct==1), x2=sum(clade.ct==2), x3=sum(clade.ct==3), x4=sum(clade.ct==4))
clades=as.matrix(clades)
#get gene counts
genes<-table(cnee.perm$accel < 0.05 & cnee.perm$tinPres > 0 & cnee.perm$tinamou > 0.1, cnee.perm$ens)[2,]
#gene level convergence (may be slow)
cnee.perm[cnee.perm$tinamou <=0.1 | cnee.perm$tinPres == 0,c(3,4,5,6)]=1
gene.conv=ddply(cnee.perm, .(ens), summarize, r.min=min(Rhea), c.min=min(Casuar), o.min=min(strCam), k.min=min(Kiwi))
gene.conv.ct=apply(gene.conv[,c("r.min", "c.min", "k.min", "o.min")], 1, function(x) sum(x < 0.05))
names(gene.conv.ct)=gene.conv$ens  
clades.real.ct=clades
genes.real.ct=genes
genes.cov.real.ct=gene.conv.ct

#permutations - version 2, this scores aceleration taking into account tinamou data first, then shuffles scores
cnee2.input = data.frame(ens=cnee.input$ens, accel=ifelse(cnee.input$accel < 0.05 & cnee.input$tinamou >= 0.1 & cnee.input$tinPres > 0, 1, 0), kiwi=ifelse(cnee.input$Kiwi < 0.05 & cnee.input$tinamou >= 0.1 & cnee.input$tinPres > 0, 1, 0), rhea=ifelse(cnee.input$Rhea < 0.05 & cnee.input$tinamou >= 0.1 & cnee.input$tinPres > 0, 1, 0), ost=ifelse(cnee.input$strCam < 0.05 & cnee.input$tinamou >= 0.1 & cnee.input$tinPres > 0, 1, 0), cas=ifelse(cnee.input$Casuar < 0.05 & cnee.input$tinamou >= 0.1 & cnee.input$tinPres > 0, 1, 0))

nreps<-1000
clade.res.2=matrix(nrow=nreps, ncol=5)
colnames(clade.res.2)=c("ct0", "ct1", "ct2", "ct3", "ct4")
gene.res.2=matrix(nrow=nreps, ncol=length(unique(cnee.input$ens)))
colnames(gene.res.2)=cnee.input$ens
genecov.res.2=gene.res.2

for (iter in 1:nreps) {
  #shuffle
  cnee.perm<-as.data.frame(apply(cnee2.input, 2, sample), stringsAsFactors=F)
  #get clade counts
  clade.ct=apply(cnee.perm[,c("rhea", "cas", "kiwi", "ost")], 1, function(x) sum(x > 0))
  clades<-data.frame(x0=sum(clade.ct==0), x1=sum(clade.ct==1), x2=sum(clade.ct==2), x3=sum(clade.ct==3), x4=sum(clade.ct==4))
  clades=as.matrix(clades)
  #get gene counts
  genes<-table(cnee.perm$accel == 1, cnee.perm$ens)[2,]
  #gene level convergence (may be slow)
  gene.conv=ddply(cnee.perm, .(ens), summarize, r.min=max(rhea), c.min=max(cas), o.min=max(ost), k.min=max(kiwi))
  gene.conv.ct=apply(gene.conv[,c("r.min", "c.min", "k.min", "o.min")], 1, function(x) sum(x > 0))
  names(gene.conv.ct)=gene.conv$ens  
  #final results
  clade.res.2[iter,]=clades
  gene.res.2[iter,]=genes
  genecov.res.2[iter,]=gene.conv.ct
}

colnames(gene.res.2)=names(genes)
colnames(genecov.res.2)=names(genes)

#real results
cnee.perm<-cnee2.input
#get clade counts
clade.ct=apply(cnee.perm[,c("rhea", "cas", "kiwi", "ost")], 1, function(x) sum(x > 0))
clades<-data.frame(x0=sum(clade.ct==0), x1=sum(clade.ct==1), x2=sum(clade.ct==2), x3=sum(clade.ct==3), x4=sum(clade.ct==4))
clades=as.matrix(clades)
#get gene counts
genes<-table(cnee.perm$accel == 1, cnee.perm$ens)[2,]
#gene level convergence (may be slow)
gene.conv=ddply(cnee.perm, .(ens), summarize, r.min=max(rhea), c.min=max(cas), o.min=max(ost), k.min=max(kiwi))
gene.conv.ct=apply(gene.conv[,c("r.min", "c.min", "k.min", "o.min")], 1, function(x) sum(x > 0))
names(gene.conv.ct)=gene.conv$ens   
clades.real.ct.2=clades
genes.real.ct.2=genes
genes.cov.real.ct.2=gene.conv.ct

#make p-values for genes
genes.sig<-data.frame(ens=rownames(sig.df), obs=sig.df$accel, total=sig.df$total, pval=1)
for (gene.to.test in genes.sig$ens) {
  obs.ct=genes.sig[genes.sig$ens==gene.to.test, "obs"]
#  pval1=sum(gene.res[,gene.to.test]>=obs.ct)/length(gene.res[,gene.to.test])
  pval2=sum(gene.res.2[,gene.to.test]>=obs.ct)/length(gene.res.2[,gene.to.test])
#  genes.sig$pval1[genes.sig$ens==gene.to.test]=pval1
  genes.sig$pval[genes.sig$ens==gene.to.test]=pval2
}

#make p-values for gene-level convergence
genes.cov.sig<-data.frame(ens=names(genes.cov.real.ct), obs=genes.cov.real.ct, pval=1)
for (gene.to.test in genes.cov.sig$ens) {
  obs.ct=genes.cov.sig[genes.cov.sig$ens==gene.to.test, "obs"]
#  pval1=sum(genecov.res[,gene.to.test]>=obs.ct)/length(genecov.res[,gene.to.test])
  pval2=sum(genecov.res.2[,gene.to.test]>=obs.ct)/length(genecov.res.2[,gene.to.test])
#  genes.cov.sig$pval1[genes.cov.sig$ens==gene.to.test]=pval1
  genes.cov.sig$pval[genes.cov.sig$ens==gene.to.test]=pval2
}

#convert to 2+, 3+, 4
twop=apply(clade.res.2[,c(3,4,5)], 1, sum)
threep=apply(clade.res.2[,c(4,5)], 1, sum)
fourp=clade.res.2[,5]

#figures for talk
hist(twop, breaks=25, ylim=c(0,100), xlim=c(0,500), col="red", las=1, main="", xlab="CNEEs with 2+ losses", cex.axis=2, cex.lab=2)
arrows(x0=sum(clades.real.ct[c(3,4,5)]),y0=20,x1=sum(clades.real.ct[c(3,4,5)]),y1=0, col="blue", lwd=6)

barplot(table(threep), ylim=c(0,700), xlim=c(0,75), col="red", las=1, main="", xlab="CNEEs with 3+ losses", cex.axis=2, cex.lab=2, xaxt="n")
axis(1, cex.axis=2)
arrows(x0=sum(clades.real.ct[c(4,5)]),y0=150,x1=sum(clades.real.ct[c(4,5)]),y1=0, col="blue", lwd=4)

barplot(table(fourp), ylim=c(0,1000), xlim=c(0,10), col="red", las=1, main="", xlab="CNEEs with 4 losses", cex.axis=2, cex.lab=2, xaxt="n")
axis(1, cex.axis=2)
arrows(x0=sum(clades.real.ct[c(5)]),y0=200,x1=sum(clades.real.ct[c(5)]),y1=0, col="blue", lwd=4)


#interesting genes
hist(genes.sig$pval[genes.sig$obs>0], breaks=seq(0,1,0.01), col="gray", freq=F, las=1, cex.axis=1.5, cex.lab=1.5, xlab="P-value", main="", ylim=c(0,6))
 abline(h=1, col="red", lty="dashed", lwd=3)

#TBX5
hist(gene.res.2[,"ENSGALG00000008253"], xlim=c(0,30), col="red", las=1, xlab="TBX accelerated CNEEs", cex.lab=1.5, cex.axis=1.5, main="")
arrows(x0=genes.sig$obs[genes.sig$ens=="ENSGALG00000008253"],y0=20,x1=genes.sig$obs[genes.sig$ens=="ENSGALG00000008253"],y1=0, col="blue", lwd=6)

#accel vs total plot
plot(x=jitter(genes.sig$obs, factor=2), y=jitter(genes.sig$total, factor=2), log="y", pch=1, col=ifelse(genes.sig$pval < 0.01, "red", "blue"), xlab="# accelerated CNEEs", ylab="Total CNEEs", las=1, cex.axis=1.5, cex.lab=1.5, cex=.9)

#add geisha expression
geisha<-read.table("~/Documents/Science/presentations/meetings/BirdMeeting 2016/geisha_expression.txt", header=T, sep="\t", quote="", comment.char="")
geisha.sub=subset(geisha, select=c("Ensembl.ID", "locations"))
geisha.sub$leg=0
geisha.sub$leg[grep("Leg", geisha.sub$location, fixed=T)]=1
geisha.sub$leg[grep("leg", geisha.sub$location, fixed=T)]=1
geisha.sub$leg[grep("wing", geisha.sub$location, fixed=T)]=1
geisha.sub$leg[grep("Wing", geisha.sub$location, fixed=T)]=1
geisha.sub=subset(geisha.sub, Ensembl.ID != "", select=c("Ensembl.ID", "leg"))
geisha.sub=geisha.sub[order(geisha.sub$leg, decreasing=T),]
geisha.clean=geisha.sub[!duplicated(geisha.sub$Ensembl.ID),]
geisha.cnee<-merge(geisha.clean, genes.sig, by.x="Ensembl.ID", by.y="ens", all.y=T)

conv.cand<-cnee.long[cnee.long$clade.ct > 2 & cnee.long$ratite.broad == 1 & cnee.long$tinPres > 0,c("id", "ens", "clade.ct", "length")]
conv.cand<-merge(conv.cand, geisha.clean, by.x="ens", by.y="Ensembl.ID", all.x=T)

##OLD CODE BELOW###
#mCE1140393
mCE1356740
mCE555774

ENSGALG00000025996
ENSGALG00000002925
ENSGALG00000016016



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
