

#OLD CODE BELOW##

#top candidates
cnee.long[cnee.long$clade.ct > 2 & cnee.long$ratite.broad == 1 & cnee.long$tinPres > 0,c("id", "ens", "clade.ct", "length")]

cnee.long$ratite.broad.strict = cnee.long$ratite.broad
cnee.long$ratite.broad.strict[cnee.long$tinPres == 0] = 0
ls

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