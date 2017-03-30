setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/hyphy_relax/")
library(data.table)

#ratites
rp<-read.table("relax_pval.all", header=F)
rk<-read.table("relax_K.all", header=F, fill=T)

rp=rp[,c(1,4,5,6,8)]
rk=rk[,c(1,4,5,6,8)]

names(rp)=c("set", "hog", "tree", "type", "pval")
names(rk)=c("set", "hog", "tree", "type", "K")
relax=merge(rp,rk)

relax=subset(relax, pval != "-------" & !is.na(K))
relax$pval=as.numeric(as.character(relax$pval))
relax=unique(relax)

#finished, good runs in relax now

#check for missing
ancrec.parsed<-fread("gunzip -c ../paml_ancrec/paml_M0_parsed.txt.gz")
ancrec.treekey<-ancrec.parsed[,c("hog", "treenum", "species_tree"), with=FALSE]
hog_info<-read.table("../all_hog_info.tsv", sep="\t", header=T)

hog_info$has_species_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree]
hog_info$has_gene_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree == F]

#merge tree info with relax
relax$tree = as.integer(sub("tree", "", relax$tree))
relax = merge(relax, ancrec.treekey, by.x=c("hog", "tree"), by.y=c("hog", "treenum"))
relax = subset(relax, species_tree)

#convert to data table
relax = as.data.table(relax)

hog_counts<-relax[,.N,by=.(hog,set)]
hogs_to_run = hog_info$hog[hog_info$has_species_tree]

write.table(hogs_to_run[!(hogs_to_run %in% hog_counts$hog[hog_counts$set=="Ratite" & hog_counts$N == 2])], file="ratite_reruns_Mar2017", quote=F, row.names = F, col.names = F)
write.table(hogs_to_run[!(hogs_to_run %in% hog_counts$hog[hog_counts$set=="VL" & hog_counts$N == 2])], file="vl_reruns_Mar2017", quote=F, row.names = F, col.names = F)
write.table(hogs_to_run[!(hogs_to_run %in% hog_counts$hog[hog_counts$set=="RND" & hog_counts$N == 2])], file="rand_reruns_Mar2017", quote=F, row.names = F, col.names = F)
write.table(hogs_to_run[!(hogs_to_run %in% hog_counts$hog[hog_counts$set=="Tinamou" & hog_counts$N == 2])], file="tinamou_reruns_Mar2017", quote=F, row.names = F, col.names = F)

## ANALYSIS

#make subset
hog_counts_all <- relax[,.N,by=hog]
hogs_to_use = hog_counts_all$hog[hog_counts_all$N == 8]

relax<-merge(relax, hog_info, by="hog")
missing_cutoff = 2
relax.clean = subset(relax, (hog %in% hogs_to_use) & dup_ct == 0 & missing_ct <= missing_cutoff, select=c("hog", "set", "type", "pval", "K"))

#add qvalues
relax.clean[,qval := p.adjust(pval, method="fdr"), by=.(set, type)]
relax.clean$sample = paste0(relax.clean$set, ".", relax.clean$type)
table(relax.clean$qval < 0.05, relax.clean$sample)

#add sig key
relax.clean$sig = as.numeric(relax.clean$qval < 0.05) * sign(1 - relax.clean$K)
#positive is increased selection, negative is relaxed selection

table(relax.clean$sig, relax.clean$sample)

tin.pos = relax.clean$hog[relax.clean$sig == 1 & relax.clean$sample == "Tinamou.tips"]
tin.neg = relax.clean$hog[relax.clean$sig == -1 & relax.clean$sample == "Tinamou.tips"]

ratite.pos = relax.clean$hog[relax.clean$sig == 1 & relax.clean$sample == "Ratite.tips"]
ratite.neg = relax.clean$hog[relax.clean$sig == -1 & relax.clean$sample == "Ratite.tips"]

relax.cmp = merge(relax.clean[sample=="VL.tips"], dn.pval, by="hog")
relax.cmp$vl.q = p.adjust(relax.cmp$vl.p, method="fdr")
fisher.test(table(relax.cmp$qval<0.01, relax.cmp$vl.q < 0.01))

#

#analysis
relax.tips<-subset(relax, species_tree==TRUE & type=="tips")
relax.tips$pval=as.numeric(as.character(relax.tips$pval))
relax.tips$qval=p.adjust(relax.tips$pval, method="fdr")

table(relax.tips$qval < 0.05, relax.tips$set)

#make sig key
relax.tips$sig = as.numeric(relax.tips$qval < 0.05) * sign(1 - relax.tips$K)
#K < 1 is relaxed, K > 1 is accelerated

table(relax.tips$set, relax.tips$sig)

#wide format
relax.tips.wide<-reshape(relax.tips[,c("hog", "set", "sig")], timevar="set", idvar="hog", direction="wide")

table(relax.tips.wide$sig.tinamou, relax.tips.wide$sig.ratite)


relax.tips$Kplot=relax.tips$K
relax.tips$Kplot[relax.tips$K < 0.01] = 0.01
relax.tips$chr=1
relax.tips$chr[relax.tips$K<0.01]=18
plot(log2(relax.tips$Kplot), -log10(relax.tips$pval), col=ifelse(relax.tips$qval<0.05,"red", "black"), xlab="log2(K value from RELAX)", ylab="-log10(p-value)", pch=relax.tips$chr)

#all
relax.all<-subset(relax, tree=="tree1" & type=="all")
relax.all$pval=as.numeric(as.character(relax.all$pval))
relax.all$qval=p.adjust(relax.all$pval, method="fdr")
relax.all$Kplot=relax.all$K
relax.all$Kplot[relax.all$K < 0.01] = 0.01
relax.all$chr=1
relax.all$chr[relax.all$K<0.01]=18
plot(log2(relax.all$Kplot), -log10(relax.all$pval), col=ifelse(relax.all$qval<0.05,"red", "black"), xlab="log2(K value from RELAX)", ylab="-log10(p-value)", pch=relax.all$chr)

#tips hogs
relaxed.hogs<-relax.tips$hog[relax.tips$K<1 & relax.tips$qval<0.05]
acc.hogs<-relax.tips$hog[relax.tips$K>1 & relax.tips$qval<0.05]

#hog->galGal
hogs<-read.table("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/homology/new_hog_list.txt")
ncbikey<-read.table("~/Downloads/CGNC_Gallus_gallus_20161020.txt", header=F, sep="\t", comment.char="", col.names=c("cgnc", "ncbi", "ensembl", "sym", "descr", "sp"), quote="")
hogs.galgal<-subset(hogs, V4=="galGal")
hogs.galgal$hog=sub("HOG2_","",hogs.galgal$V1,fixed=T)
hogs.galgal<-merge(hogs.galgal, ncbikey, all.x=T, all.y=F, by.x="V3", by.y="ncbi")

relaxed.hogs.galgal<-subset(hogs.galgal, hog %in% relaxed.hogs)
write.table(relaxed.hogs.galgal$sym, file="relaxed.galgal.sym", quote=F, row.names=F, col.names=F)
write.table(relaxed.hogs.galgal$ensembl, file="relaxed.galgal.ens", quote=F, row.names=F, col.names=F)
write.table(relaxed.hogs.galgal$V3, file="relaxed.galgal.ncbi", quote=F, row.names=F, col.names=F)

accel.hogs.galgal<-subset(hogs.galgal, hog %in% acc.hogs)
write.table(accel.hogs.galgal$V3, file="acc.galgal.ncbi", quote=F, row.names=F, col.names=F)
