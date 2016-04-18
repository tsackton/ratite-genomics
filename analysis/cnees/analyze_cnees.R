library(data.table)

#load required inputs

ce_all.annot<-fread("/Volumes/LaCie/Projects/Current/ratites/final/ce_final/final_ces.tree2.bed.final_annotation.tsv", header=T)
ce_nor.annot<-fread("/Volumes/LaCie/Projects/Current/ratites/final/ce_final/final_ces_noratite.tree2.bed.final_annotation.tsv", header=T)
load("/Volumes/LaCie/Projects/Current/ratites/final/accelTests/all.accel.processed.Rdata")
alloc.col(all.accel)
#this is not ideal
load("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/cnee_presence/cnee.pres.full.Rlist")
cnee.pres.all<-cnee.pres
load("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/cnee_presence/cnee.pres.rr_only.Rlist")
cnee.pres.rronly<-cnee.pres
rm(cnee.pres)

#now have everything we need for filtering

#first step -- make subsets of annotation files just for >=50 elements
ce.annot<-rbind(ce_nor.annot[id %in% all.accel[runtype=="rrOnlyCNEEs",name],], ce_all.annot[length>50 & class !="CDS" & class != "exonic_non_CDS",])
#ce.annot has all the CNEEs we care about

#remove other annotation objects
rm(ce_all.annot)
rm(ce_nor.annot)

#second step, process presence/absence data
#convert to data.table for processing efficiency
all.pres<-do.call("rbind", cnee.pres.all)
all.pres.dt<-as.data.table(all.pres)
nr.pres<-do.call("rbind", cnee.pres.rronly)
nr.pres.dt<-as.data.table(nr.pres)
cnee.pres<-rbind(all.pres.dt, nr.pres.dt)
cnee.pres<-cnee.pres[name %in% ce.annot$id,]
cnee.pres<-merge(cnee.pres, ce.annot[,c("id", "length"), with=FALSE], by.x="name", by.y="id")
cnee.pres$rel.len = cnee.pres$target.length/cnee.pres$length

##SET CUTOFFS##
length.cutoff<-0.5
dup.cutoff<-1
#############

cnee.pres$present=ifelse(cnee.pres$rel.len >= length.cutoff, T, F)
cnee.pres$duplicated=ifelse(cnee.pres$count > dup.cutoff, T, F)

#now make the final presence file
cnee.pres.final<-as.data.frame(table(cnee.pres$duplicated,cnee.pres$name), stringsAsFactors=F)
cnee.pres.final=cnee.pres.final[cnee.pres.final$Var1=="TRUE", c("Var2", "Freq")]
names(cnee.pres.final)=c("name", "dupcount")

species.list<-c("allMis", "aptFor", "balReg", "chaVoc", "corBra", "droNov", "lepDis", "nipNip", "pygAde", "taeGut", "allSin", "aptHaa", "calAnn", "cheMyd", "croPor", "eudEle", "fulGla", "melGal", "notPer", "rheAme", "tinGut", "anaPla", "aptOwe", "casCas", "chrPic", "cryCin", "falPer", "gavGan", "melUnd", "picPub", "rhePen", "anoCar", "aptRow", "chaPel", "colLiv", "cucCan", "ficAlb", "halLeu", "mesUni", "pseHum", "strCam")

#add presence to cnee.data
for (species in species.list) {
  newcol=paste0(species, "Pres")
  cnee.pres.final[,newcol]=cnee.pres.final$name %in% cnee.pres[sp==species & present==TRUE, name]
}

cnee.pres.final$tinPres = apply(cnee.pres.final[,c("tinGutPres", "eudElePres", "notPerPres", "cryCinPres")], 1, sum)
cnee.pres.final$ratitePres = apply(cnee.pres.final[,c("aptHaaPres", "aptRowPres", "aptOwePres", "strCamPres", "rhePenPres", "rheAmePres", "casCasPres", "droNovPres")], 1, sum)
cnee.pres.final$nonAvesPres = apply(cnee.pres.final[,c("anoCarPres", "chrPicPres", "cheMydPres", "allSinPres", "allMisPres", "gavGanPres", "croPorPres")], 1, sum)
cnee.pres.final$AvesPres = apply(cnee.pres.final[,c("aptForPres", "balRegPres", "chaVocPres", "corBraPres", "droNovPres", "lepDisPres", "nipNipPres", "pygAdePres", "taeGutPres", "aptHaaPres", "calAnnPres",  "eudElePres", "fulGlaPres", "melGalPres", "notPerPres", "rheAmePres", "tinGutPres", "anaPlaPres", "aptOwePres", "casCasPres", "cryCinPres", "falPerPres", "melUndPres", "picPubPres", "rhePenPres","aptRowPres", "chaPelPres", "colLivPres", "cucCanPres", "ficAlbPres", "halLeuPres", "mesUniPres", "pseHumPres", "strCamPres")], 1, sum)+1

#remove unneeded cnee presence objects
rm(newcol)
rm(nr.pres)
rm(species)
rm(species.list)
rm(all.pres.dt)
rm(all.pres)
rm(cnee.pres.all)
rm(cnee.pres.rronly)
rm(nr.pres.dt)

#now clean up acceleration data to limit to the proper CNEEs
accel.final<-all.accel[name %in% ce.annot$id,]
accel.final$align="with_moa"
accel.final$align[accel.final$runtype=="rrOnlyCNEEs" | accel.final$runtype=="orig"]="wga"

accel.final$result.type=paste(accel.final$align, accel.final$model, sep=".")
#remove runtype, align, and model columns, just keep result.type
accel.final[,c("runtype","align","model") := NULL]

#add variables
#subtree rate is relative branch length of subtree (defined by group column)
#twop is p-value converted to a two-tailed version
accel.final$subtree.rate = accel.final$alt_subscale * accel.final$alt_scale
accel.final$twop=accel.final$pval*2
accel.final$twop[accel.final$twop==0]=1e-5
accel.final$twop[accel.final$twop>1]=1
accel.final[,qval.raw:=p.adjust(pval, method="fdr", n=285205), by=c("result.type", "group")]

#filter to remove CNEEs with too much missing or duplicated data
cnees.good<-cnee.pres.final[cnee.pres.final$dupcount <= 3 & cnee.pres.final$AvesPres >= 27 & cnee.pres.final$tinPres >= 3, c("name")]

##ADD Q-VALUES TO FILTERED DATASET##
accel.filtered<-accel.final[name %in% cnees.good,]

#testing some fancy-pants data.table syntax for update column by group
accel.filtered[,qval1:=p.adjust(pval, method="fdr", n=270414), by=c("result.type", "group")]
accel.filtered[,qval2:=p.adjust(twop, method="fdr", n=270414), by=c("result.type", "group")]

#make a rate matrix and a qvalue matrix
accel.qval<-dcast(accel.filtered, name ~ group + result.type, value.var="qval1")
accel.rate<-dcast(accel.filtered, name ~ group + result.type, value.var="subtree.rate")

##SENSITIVITY ANALYSIS##
##need to move from labmeeting.R##

#generate data for evo-devo work
accel.clade<-dcast(accel.filtered, name + result.type ~ group, value.var="qval1")
accel.clade$filtQ=apply(accel.clade[,c("BasalPaleo","BasalPaleo_noR","Tinamou"),with=F], 1, min, na.rm=T)
accel.clade$cas.acc=ifelse(accel.clade$Casuar <= 0.05 & accel.clade$filtQ >= 0.25, 1, 0)
accel.clade$kiwi.acc=ifelse(accel.clade$Kiwi <= 0.05 & accel.clade$filtQ >= 0.25, 1, 0)
accel.clade$moa.acc=ifelse(!is.na(accel.clade$anoDid) & accel.clade$anoDid <= 0.05 & accel.clade$filtQ >= 0.25, 1, 0)
accel.clade$ostrich.acc=ifelse(accel.clade$strCam <= 0.05 & accel.clade$filtQ >= 0.25, 1, 0)
accel.clade$rhea.acc=ifelse(accel.clade$Rhea <= 0.05 & accel.clade$filtQ >= 0.25, 1, 0)

#get clade count
accel.clade[,clade.ct:=rhea.acc+cas.acc+kiwi.acc+moa.acc+ostrich.acc]
              
#make into wide format
accel.wide<-dcast(accel.clade, name ~ result.type, value.var="clade.ct")
accel.wide$max=apply(accel.wide[,c("with_moa.neut_ver3","with_moa.neut_ver2","with_moa.neut_ver1","wga.neut_ver2"),with=F], 1, max, na.rm=T)
accel.wide$min=apply(accel.wide[,c("with_moa.neut_ver3","with_moa.neut_ver2","with_moa.neut_ver1","wga.neut_ver2"),with=F], 1, min, na.rm=T)

accel.wide<-merge(accel.wide, ce_all[,c("id", "best_ens", "length", "class"),with=F], by.x="name", by.y="id", all.x=T, all.y=F)
accel.clade<-merge(accel.clade, ce_all[,c("id", "best_ens", "length", "class"),with=F], by.x="name", by.y="id", all.x=T, all.y=F)

#write out
write.table(accel.wide, file="clade_count_updated.tsv", sep="\t", quote=F, row.names=F)
write.table(accel.clade, file="clade_count_updated_raw.tsv", sep="\t", quote=F, row.names=F)

#presence/absence
#get all available acceleration tests for palaeognath clade as a whole
accel.pal<-dcast(accel.final[group == "Tinamou" | group=="BasalPaleo" | group=="BasalPaleo_noR",], name ~ result.type + group, value.var="qval.raw")
accel.pal$filtQ = apply(accel.pal[,2:10,with=F], 1, min, na.rm=T)

cnee.loss<-merge(cnee.pres.final, accel.pal[,c("name", "filtQ"),with=F],by="name", all.x=T, all.y=F)
cnee.loss.hc<-subset(cnee.loss, filtQ>=0.25 & tinPres == 4)
cnee.loss.hc$nonRatiteLoss = 35 - cnee.loss.hc$AvesPres - (8 - cnee.loss.hc$ratitePres)

cnee.loss.hc<-merge(cnee.loss.hc, ce_all[,c("id", "best_ens", "length", "class"),with=F], by.x="name", by.y="id", all.x=T, all.y=F)
write.table(cnee.loss.hc[cnee.loss.hc$nonRatiteLoss==0 & cnee.loss.hc$ratitePres < 8, c("name", "best_ens", "length", "class", "ratitePres", "nonAvesPres", "strCamPres", "rheAmePres", "rhePenPres", "aptHaaPres", "aptOwePres", "aptRowPres", "droNovPres", "casCasPres")], file="ratite_cnee_losses.tsv", sep="\t", row.names=F, quote=F)
