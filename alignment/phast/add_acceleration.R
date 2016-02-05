#read in ce annotation
ce.annot <- read.table("./final_beds/ce_annotation.tsv", header=T, stringsAsFactors = F)

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

#length subset
cnee.long = subset(cnee, length > 50)

#write out ids
write.table(cnee.long[,c("id")], "cnee_long.ids", row.names=F, quote=F, sep="\t", col.names=Ftabl)

#get elements present in each species

psp <- c("allMis", "aptFor", "balReg", "chaVoc", "corBra", "droNov", "lepDis", "nipNip", "pygAde", "taeGut", "allSin", "aptHaa", "calAnn", "cheMyd", "croPor", "eudEle", "fulGla", "melGal", "notPer", "rheAme", "tinGut", "anaPla", "aptOwe", "casCas", "chrPic", "cryCin", "falPer", "gavGan", "melUnd", "picPub", "rhePen", "anoCar", "aptRow", "chaPel", "colLiv", "cucCan", "ficAlb", "halLeu", "mesUni", "pseHum", "strCam")
cnee.pres<-list()
for (group in psp) {
  cnee.pres[[group]]<-read.table(paste0("final_accel/", group, ".bed", sep=""), header=F, sep="\t", comment.char="")
  cnee.pres[[group]]$sp = group
}

#add presence to cnee.long
for (sp in psp) {
  newcol=paste0(sp, "Pres")
  cnee.long[,newcol]=cnee.long$id %in% cnee.pres[[sp]]$V4
}

cnee.long$tinPres = apply(cnee.long[,c("tinGutPres", "eudElePres", "notPerPres", "cryCinPres")], 1, sum)
cnee.long$ratitePres = apply(cnee.long[,c("aptHaaPres", "aptRowPres", "aptOwePres", "strCamPres", "rhePenPres", "rheAmePres", "casCasPres", "droNovPres")], 1, sum)
cnee.long$nonAvesPres = apply(cnee.long[,c("anoCarPres", "chrPicPres", "cheMydPres", "allSinPres", "allMisPres", "gavGanPres", "croPorPres")], 1, sum)
cnee.long$AvesPres = apply(cnee.long[,c("aptForPres", "balRegPres", "chaVocPres", "corBraPres", "droNovPres", "lepDisPres", "nipNipPres", "pygAdePres", "taeGutPres", "aptHaaPres", "calAnnPres",  "eudElePres", "fulGlaPres", "melGalPres", "notPerPres", "rheAmePres", "tinGutPres", "anaPlaPres", "aptOwePres", "casCasPres", "cryCinPres", "falPerPres", "melUndPres", "picPubPres", "rhePenPres","aptRowPres", "chaPelPres", "colLivPres", "cucCanPres", "ficAlbPres", "halLeuPres", "mesUniPres", "pseHumPres", "strCamPres")], 1, sum)+1

#write out cnee analysis data
write.table(cnee.long, "cnee_long_info.tsv", sep="\t", quote=F, row.names=F)
