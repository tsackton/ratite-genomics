indir<-"/Volumes/LaCie/Projects/Current/ratites/final/accelTests/withMoa_neognath"
setwd(indir)

chaVoc<-read.table("chaVoc.out_neut_ver3.results", header=T, stringsAsFactors=F)
colLiv<-read.table("colLiv.out_neut_ver3.results", header=T, stringsAsFactors=F)
halLeu<-read.table("halLeu.out_neut_ver3.results", header=T, stringsAsFactors=F)
nipNip<-read.table("nipNip.out_neut_ver3.results", header=T, stringsAsFactors=F)
taeGut<-read.table("taeGut.out_neut_ver3.results", header=T, stringsAsFactors=F)

neoaccl<-data.frame(ce=chaVoc$name, chaVoc=p.adjust(chaVoc$pval, method="fdr"), colLiv=p.adjust(colLiv$pval, method="fdr"), halLeu=p.adjust(halLeu$pval, method="fdr"), nipNip=p.adjust(nipNip$pval, method="fdr"), taeGut=p.adjust(taeGut$pval, method="fdr"), stringsAsFactors=F)

neoaccl.good = subset(neoaccl, neoaccl$ce %in% cnees.good, select=c(2,3,4,5,6))
neoaccl.good$clade.ct = apply(neoaccl.good, 1, function(x) sum(x < 0.05))

setwd("~/Projects/")

