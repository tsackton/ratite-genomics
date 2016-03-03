#slow and clunky R code to read psl file and parse it to estimate coverage and duplications
library(plyr)
library(reshape)

psp <- c("allMis", "aptFor", "balReg", "chaVoc", "corBra", "droNov", "lepDis", "nipNip", "pygAde", "taeGut", "allSin", "aptHaa", "calAnn", "cheMyd", "croPor", "eudEle", "fulGla", "melGal", "notPer", "rheAme", "tinGut", "anaPla", "aptOwe", "casCas", "chrPic", "cryCin", "falPer", "gavGan", "melUnd", "picPub", "rhePen", "anoCar", "aptRow", "chaPel", "colLiv", "cucCan", "ficAlb", "halLeu", "mesUni", "pseHum", "strCam")
cnee.pres<-list()
for (group in psp) {
  psl<-read.table(paste0("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/alignment/cnee_presence/", group, ".psl", sep=""), header=F, sep="\t", comment.char="", stringsAsFactors=F)
  psl$q.length=psl$V14-psl$V13
  psl$t.length=apply(psl[,"V20", drop=F], 1, function(x) sum(as.numeric(unlist(strsplit(x, ",", fixed=T)))))  
  dupct<-as.data.frame(table(psl$V1), stringsAsFactors=F)
  tlen<-data.frame(Var1=psl$V1, len=psl$t.length, stringsAsFactors=F)
  output<-merge(dupct, tlen, all=T, by="Var1")
  names(output)=c("name", "count", "target.length")
  cnee.pres[[group]]=ddply(output, .(name), summarize, count=max(count), target.length=sum(target.length))
  cnee.pres[[group]]$sp=group 
}

#write out list as a backup
save(cnee.pres, file="cnee.pres.Rlist")

#convert to data.table for processing efficiency
all.pres<-do.call("rbind", cnee.pres)
library(data.table)
all.pres.dt<-as.data.table(all.pres)
rm(all.pres)

#okay now we have a nice data table with the target length and copy number
#need to add the query length to each row

ce.len<-as.data.table(ce.annot[,c("id", "length")])

all.pres.dt<-merge(all.pres.dt, ce.len, by.x="name", by.y="id", all.x=T, all.y=F)
all.pres.dt$rel.len = all.pres.dt$target.length/all.pres.dt$length

#cutoff sets the required coverage to call as present
cutoff<-0.5
#dupmax sets the tolerated number of duplications prior to calling as duplicated
dupmax<-1

all.pres.dt$present=ifelse(all.pres.dt$rel.len >= cutoff, T, F)
all.pres.dt$duplicated=ifelse(all.pres.dt$count > dupmax, T, F)

dup.ct.percnee<-as.data.frame(table(all.pres.dt[duplicated==T,name]))

#now we do the final calculations
cnee.data<-ce.annot[,c("id", "length", "in.intersection", "class", "cnee", "best_ncbi", "best_ens")]

cnee.data<-merge(cnee.data, dup.ct.percnee, by.x="id", by.y="Var1", all.x=T)
cnee.data$Freq[is.na(cnee.data$Freq)]=0

#add presence to cnee.data
for (species in psp) {
  newcol=paste0(species, "Pres")
  cnee.data[,newcol]=cnee.data$id %in% all.pres.dt[sp==species & present==TRUE, name]
}


cnee.data$tinPres = apply(cnee.data[,c("tinGutPres", "eudElePres", "notPerPres", "cryCinPres")], 1, sum)
cnee.data$ratitePres = apply(cnee.data[,c("aptHaaPres", "aptRowPres", "aptOwePres", "strCamPres", "rhePenPres", "rheAmePres", "casCasPres", "droNovPres")], 1, sum)
cnee.data$nonAvesPres = apply(cnee.data[,c("anoCarPres", "chrPicPres", "cheMydPres", "allSinPres", "allMisPres", "gavGanPres", "croPorPres")], 1, sum)
cnee.data$AvesPres = apply(cnee.data[,c("aptForPres", "balRegPres", "chaVocPres", "corBraPres", "droNovPres", "lepDisPres", "nipNipPres", "pygAdePres", "taeGutPres", "aptHaaPres", "calAnnPres",  "eudElePres", "fulGlaPres", "melGalPres", "notPerPres", "rheAmePres", "tinGutPres", "anaPlaPres", "aptOwePres", "casCasPres", "cryCinPres", "falPerPres", "melUndPres", "picPubPres", "rhePenPres","aptRowPres", "chaPelPres", "colLivPres", "cucCanPres", "ficAlbPres", "halLeuPres", "mesUniPres", "pseHumPres", "strCamPres")], 1, sum)+1


