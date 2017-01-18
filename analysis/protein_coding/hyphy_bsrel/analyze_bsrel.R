setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/hyphy_bsrel/")
library(data.table)

#make hog<->galGal key
hogs<-read.table("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/homology/new_hog_list.txt")
ncbikey<-read.table("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/annotation/galGalAnnot/CGNC_Gallus_gallus_20161020.txt", header=F, sep="\t", comment.char="", col.names=c("cgnc", "ncbi", "ensembl", "sym", "descr", "sp"), quote="")
hogs.galgal<-subset(hogs, V4=="galGal")
hogs.galgal$hog=sub("HOG2_","",hogs.galgal$V1,fixed=T)
hogs.galgal<-merge(hogs.galgal, ncbikey, all.x=T, all.y=F, by.x="V3", by.y="ncbi")

merged = read.table("merged_res.txt", fill=T)
names(merged)=c("class", "tree", "hog", "tsel", "nsel", "tnon", "nnon")

#load tree key from paml_ancrec 
ancrec.parsed<-fread("gunzip -c ../paml_ancrec/ancrec_parsed.out.gz")
paml.treekey<-ancrec.parsed[,c("hog", "treenum", "species_tree"), with=FALSE]
paml.treekey$tree = paste0("tree", paml.treekey$treenum)
merged=merge(merged, paml.treekey, by.x=c("hog", "tree"), by.y=c("hog", "tree"))
merged = subset(merged, species_tree==TRUE)

#note that this is just the ~8700 hogs with no duplication and just the species tree runs

merged$total_sel = merged$tsel + merged$nsel
merged$target_prop = merged$tsel / merged$total_sel
table(merged$total_sel, merged$class)

#add annotation
bsrel<-merge(merged, hogs.galgal, by="hog", all.x=T, all.y=F)
bsrel$tot_lin = bsrel$tsel + bsrel$nsel + bsrel$tnon + bsrel$nnon
bsrel$target_lin = bsrel$tsel + bsrel$tnon

bsrel.clean = subset(bsrel, target_lin > 0)

plot(sort(bsrel.clean$target_prop[bsrel.clean$class=="ratites"]), type='l', col="red", ylab="Proportion of Selected Branches in Target Class")
lines(sort(bsrel.clean$target_prop[bsrel.clean$class=="vl"]), type="l", col="darkgreen")
lines(sort(bsrel.clean$target_prop[bsrel.clean$class=="rand1"]), type="l", col="gray")
lines(sort(bsrel.clean$target_prop[bsrel.clean$class=="rand2"]), type="l", col="gray")
legend("topleft", legend=c("'Random'", "Ratites", "Vocal Learners"), col=c("gray", "red",  "darkgreen"), lwd=2)

table(bsrel$tsel, bsrel$target_prop == 1, bsrel$class, useNA="ifany")

barplot(table(bsrel$tsel[bsrel$nsel == 0 & bsrel$tsel > 1], bsrel$class[bsrel$nsel == 0 & bsrel$tsel > 1]), ylab="Number of Genes Uniquely Selected in >1 Target Lineages")
legend("topleft", legend=c("2", "3"), col=c("gray20", "gray80"), pch=15, cex=2)


#ratite specific selection
table(bsrel$target_prop==1, bsrel$tot_lin, bsrel$class)
table(bsrel$target_prop==1, bsrel$tsel>1, bsrel$class)

#gene symbol list
write.table(bsrel$sym[bsrel$class=="ratites" & bsrel$target_prop==1], file="ratite.selected", quote=F, row.names=F, col.names=F)
write.table(bsrel$sym[bsrel$class=="ratites" & bsrel$target_prop<1 & bsrel$total_sel > 0], file="ratite.background", quote=F, row.names=F, col.names=F)
write.table(bsrel$sym[bsrel$class=="vl" & bsrel$target_prop==1], file="vl.selected", quote=F, row.names=F, col.names=F)
write.table(bsrel$sym[bsrel$class=="vl" & bsrel$target_prop<1 & bsrel$total_sel > 0], file="vl.background", quote=F, row.names=F, col.names=F)

