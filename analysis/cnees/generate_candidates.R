#generate candidates 
#currently requires the analyze_cnees.R script to be run first
#will refactor to load data eventually

library(ape)

#make into wide format
accel.wide<-dcast(accel.clade, name ~ result.type, value.var=c("clade.ct", "clade.ct.2"))
accel.wide$max.1=apply(accel.wide[,c("clade.ct_with_moa.neut_ver3","clade.ct_with_moa.neut_ver2","clade.ct_with_moa.neut_ver1","clade.ct_wga.neut_ver1", "clade.ct_wga.neut_ver2", "clade.ct_wga.neut_ver3"),with=F], 1, max, na.rm=T)
accel.wide$min.1=apply(accel.wide[,c("clade.ct_with_moa.neut_ver3","clade.ct_with_moa.neut_ver2","clade.ct_with_moa.neut_ver1","clade.ct_wga.neut_ver1", "clade.ct_wga.neut_ver2", "clade.ct_wga.neut_ver3"),with=F], 1, min, na.rm=T)
accel.wide$max.2=apply(accel.wide[,c("clade.ct.2_with_moa.neut_ver3","clade.ct.2_with_moa.neut_ver2","clade.ct.2_with_moa.neut_ver1","clade.ct.2_wga.neut_ver1", "clade.ct.2_wga.neut_ver2", "clade.ct.2_wga.neut_ver3"),with=F], 1, max, na.rm=T)
accel.wide$min.2=apply(accel.wide[,c("clade.ct.2_with_moa.neut_ver3","clade.ct.2_with_moa.neut_ver2","clade.ct.2_with_moa.neut_ver1","clade.ct.2_wga.neut_ver1", "clade.ct.2_wga.neut_ver2", "clade.ct.2_wga.neut_ver3"),with=F], 1, min, na.rm=T)
accel.wide<-merge(accel.wide, ce.annot[,c("id", "best_ens", "length", "class"),with=F], by.x="name", by.y="id", all.x=T, all.y=F)

#add score from cons score analysis
accel.wide<-merge(accel.wide, clust, by.x="name", by.y="elementID", all.x=T, all.y=F)

#get minimum acceleration q-value from accel.clade
accel.pal<-dcast(accel.final[group == "Tinamou" | group=="BasalPaleo" | group=="BasalPaleo_noR",], name ~ result.type + group, value.var="qval.raw")
accel.pal$filtQ = apply(accel.pal[,2:10,with=F], 1, min, na.rm=T)
accel.wide<-merge(accel.wide, as.data.frame(accel.pal)[,c("name","filtQ")], all.x=T, all.y=F, by="name")
accel.wide<-merge(accel.wide, cnee.pres.final, all.x=T, all.y=F, by="name")
accel.wide$neognathLoss = 21 - (accel.wide$AvesPres - accel.wide$ratitePres - accel.wide$tinPres - accel.wide$aptForPres - accel.wide$pygAdePres)
#ratite candidate set generation
#define best candidate set as:
#1 (clustType == convFlightless | clustType == convFlightless) & (min.2 >= 2) 

ratite.candidates=subset(accel.wide, (clustType == "convFlightless" | clustType == "convRatite" | is.na(clustType)) & (min.2 >= 2) & (filtQ > 0.25) & (neognathLoss == 0) & (tinPres == 4))

#add symbols
ens.info<-read.table("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/annotation/galGalAnnot//galGal_ens_sym_04202016.txt", header=T, sep="\t", stringsAsFactors=F, quote="")

ratite.candidates<-merge(ratite.candidates, ens.info, by.x="best_ens", by.y="Ensembl.Gene.ID", all.x=T, all.y=F)

write.table(ratite.candidates[,c("name", "best_ens", "Associated.Gene.Name", "Description", "Chromosome.Name", "max.1", "max.2", "min.1", "min.2", "length", "clustType", "ratitePres", "nonAvesPres"), with=F], file="ratite_candidates_final.tsv", sep="\t", quote=F, row.names=F)

#read in and plot all trees associated with candidates
load("/Volumes/LaCie/Projects/Current/ratites/final/cnee_branchlen/alltrees.Rlist")
cand.trees<-data.frame(path=alltrees, name=character(284001), stringsAsFactors=F)
cand.trees$name = sub(".*RAxML_result\\.", "", cand.trees$path, perl=T)
cand.trees$name<-sub("_ver\\d+", "", cand.trees$name, perl=T)
cand.trees=subset(cand.trees, name %in% ratite.candidates$name)

#get index
fasta.index<-fread("/Volumes/LaCie/Projects/Current/ratites/final/cnee_branchlen/fasta_records_list", header=F, sep="\t")
fasta.index$V1=sub("input_fastas_allspecies_cnees_aligned_no-galGal-gaps/batch\\d+\\/", "", fasta.index$V1, perl=T)
fasta.index$V1=sub("\\.fasta:", "", fasta.index$V1, perl=T)

pdf(file="candidate_cnee_trees.pdf", width=8, height=12)
for (i in 1:length(cand.trees$name)) {
  intree<-read.tree(cand.trees$path[i])
  treename<-cand.trees$name[i]
  ttd<-intree$tip.label[!(intree$tip.label %in% fasta.index$V2[fasta.index$V1==cand.trees$name[i]])]
  plot(drop.tip(root(intree, outgroup="anoCar"), tip=unique(c(ttd, "anoCar", "cheMyd", "chrPic", "gavGan", "croPor", "allSin", "allMis"))), no.margin=T, cex=0.85)
  title(main=treename, line=-2)
}
dev.off()
