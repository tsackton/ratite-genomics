#slow and clunky R code to read psl file and parse it to estimate coverage and duplications
library(plyr)
library(reshape)
path="/Volumes/LaCie/Projects/Current/ratites/final/cnee_presence/rr_only/"
outfile="cnee.pres.rr_only.Rlist"

psp <- c("allMis", "aptFor", "balReg", "chaVoc", "corBra", "droNov", "lepDis", "nipNip", "pygAde", "taeGut", "allSin", "aptHaa", "calAnn", "cheMyd", "croPor", "eudEle", "fulGla", "melGal", "notPer", "rheAme", "tinGut", "anaPla", "aptOwe", "casCas", "chrPic", "cryCin", "falPer", "gavGan", "melUnd", "picPub", "rhePen", "anoCar", "aptRow", "chaPel", "colLiv", "cucCan", "ficAlb", "halLeu", "mesUni", "pseHum", "strCam")
cnee.pres<-list()
for (group in psp) {
  psl<-read.table(paste0(path, group, ".psl", sep=""), header=F, sep="\t", comment.char="", stringsAsFactors=F)
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
saveRDS(cnee.pres, file=outfile)

