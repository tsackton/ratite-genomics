load("~/Downloads/models.Rdata")
load("~/Downloads/modelscores.Rdata")

birds <- c("mesUni","pygAde","aptFor","nipNip","aptHaa","aptOwe","aptRow","casCas","droNov","rheAme","rhePen","strCam", "galGal","melGal","taeGut","melUnd","ficAlb","pseHum","colLiv","falPer","anaPla","fulGla","lepDis","corBra","picPub","calAnn","chaVoc","cucCan","balReg","halLeu","chaPel","cryCin","tinGut","eudEle","notPer")
ratites    <-c("aptHaa","aptOwe","aptRow","casCas","droNov","rheAme","rhePen","strCam")  # flightless only
flightless <- c("mesUni","pygAde","aptFor","nipNip","aptHaa","aptOwe","aptRow","casCas","droNov","rheAme","rhePen","strCam")
flight <- c("galGal","melGal","taeGut","melUnd","ficAlb","pseHum","colLiv","falPer","anaPla","fulGla","lepDis","corBra","picPub","calAnn","chaVoc","cucCan","balReg","halLeu","chaPel","cryCin","tinGut","eudEle","notPer")
ancestral <- c("croPor", "gavGan", "anoCar","chrPic","allMis","cheMyd","allSin")

vocal <- c('calAnn', 'corBra', 'serCan', 'geoFor', 'melUnd', 'pseHum', 'taeGut', 'ficAlb')

sdata <- data.frame(taxa=c(ratites, birds[is.na(match(birds, ratites))]), is.flightless=c(rep(1, length(ratites)), rep(0, length(birds[is.na(match(birds, ratites))]))))

sdata <- cbind.data.frame(sdata, modelscores[!is.na(match(row.names(modelscores), sdata[,1])),])

sdata$is.vl <- as.numeric(sdata$taxa %in% vocal)
pvals<-numeric(length(modelset))
coef<-numeric(length(modelset))

for(j in 1:length(modelset)) {
  i<-modelset[j]
  test <- with(sdata, wilcox.test(get(i) ~ is.flightless))
  pvals[j] <- test$p.value
}


pvals.vl<-numeric(length(modelset))
coef.vl<-numeric(length(modelset))

for(j in 1:length(modelset)) {
  i<-modelset[j]
  test <- with(sdata, wilcox.test(get(i) ~ is.vl))
  pvals.vl[j] <- test$p.value
}


cor.test(sdata[3,3:11178],sdata[4,3:11178], method="sp")
cor.test(unlist(sdata[22,3:11278]),unlist(sdata[25,3:11278]), method="sp")


#try standarizing

scaledata<-scale(modelscores, center=T, scale=T)
cor(t(scaledata), use="na.or.complete")

rownames(scaledata)
scaledf<-as.data.frame(scaledata)
scaledf$is.flightless<-as.numeric(rownames(scaledf) %in% flightless)

pvals.rs<-numeric(length(modelset))
for(j in 1:length(modelset)) {
  i<-modelset[j]
  pval<-NA
  try(pval <- with(scaledf, wilcox.test(get(i) ~ is.flightless))$p.value)
  pvals.rs[j] <- pval
}


scaledf$is.vl<-as.numeric(rownames(scaledf) %in% vocal)
pvals.vs<-numeric(length(modelset))
for(j in 1:length(modelset)) {
  i<-modelset[j]
  pval<-NA
  try(pval <- with(scaledf, wilcox.test(get(i) ~ is.vl))$p.value)
  pvals.vs[j] <- pval
}

str(scaledf)
boxplot(unlist(scaledf[1,1:10000]),unlist(scaledf[2,1:10000]),unlist(scaledf[3,1:10000]),unlist(scaledf[4,1:10000]),unlist(scaledf[5,1:10000]),unlist(scaledf[6,1:10000]),unlist(scaledf[7,1:10000]),unlist(scaledf[8,1:10000]),unlist(scaledf[9,1:10000]),unlist(scaledf[10,1:10000]),unlist(scaledf[11,1:10000]),unlist(scaledf[12,1:10000]),outline=F)

scaledf.t<-t(scaledf)
scaledf.tscale<-scale(scaledf.t)

boxplot(scaledf.tscale[,1:42])

finaldf<-as.data.frame(t(scaledf.tscale))
finaldf$is.flightless<-as.numeric(rownames(finaldf) %in% flightless)
finaldf$is.vl<-as.numeric(rownames(finaldf) %in% vocal)

pvals.vss<-numeric(length(modelset))
for(j in 1:length(modelset)) {
  i<-modelset[j]
  pval<-NA
  try(pval <- with(finaldf, wilcox.test(get(i) ~ is.vl))$p.value)
  pvals.vss[j] <- pval
}


pvals.rss<-numeric(length(modelset))
for(j in 1:length(modelset)) {
  i<-modelset[j]
  pval<-NA
  try(pval <- with(finaldf, wilcox.test(get(i) ~ is.flightless))$p.value)
  pvals.rss[j] <- pval
}
