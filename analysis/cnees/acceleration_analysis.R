setwd("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/")

#read moa raw data
accel.classes<-c("Casuar", "Rhea", "Kiwi", "tinamou", "allRatite", "aptOwe", "aptHaa", "aptRow", "casCas", "droNov", "rheAme", "rhePen", "strCam", "anoDid", "basalPaleo")
neut.mods<-c("neut_ver1", "neut_ver2", "neut_ver3")
accel.res<-list()
for (neut in neut.mods) {
  for (group in accel.classes) {
    accel.res[[neut]][[group]]<-read.table(paste0("moa/", group, ".out_", neut, ".results", sep=""), header=T, sep="\t", comment.char="")
    accel.res[[neut]][[group]]$qval = p.adjust(accel.res[[neut]][[group]]$pval, method="fdr")
    accel.res[[neut]][[group]]$group = group
    accel.res[[neut]][[group]]$model = neut    
  }
}

moa.final<-data.frame(name=character(0), null_scale=numeric(0), alt_scale=numeric(0), alt_subscale=numeric(0), lnlratio=numeric(0), pval=numeric(0), qval=numeric(0), group=character(0), model=character(0))

for (neut in neut.mods) {
  accel.temp<-accel.res[[neut]]
  accel.temp2<-do.call("rbind", accel.temp)
  moa.final<-rbind(moa.final, accel.temp2)
}

#read original data

#load info
accel.classes<-c("Casuar", "Rhea", "Kiwi", "tinamou", "ratite", "aptOwe", "aptHaa", "aptRow", "casCas", "droNov", "rheAme", "rhePen", "strCam")
accel.res<-list()
for (group in accel.classes) {
  accel.res[[group]]<-read.table(paste0("phast/final_accel/all_", group, ".out", sep=""), header=T, sep="\t", comment.char="")
}

for (group in accel.classes) {
  accel.res[[group]]$qval = p.adjust(accel.res[[group]]$pval, method="fdr")
  accel.res[[group]]$group = group
}

accel.orig<-do.call("rbind", accel.res)
accel.new<-moa.final

accel.orig$subtree_rate = accel.orig$alt_scale * accel.orig$alt_subscale
accel.new$subtree_rate = accel.new$alt_scale * accel.new$alt_subscale

accel.orig.wide<-spread(subset(accel.orig, select=c("name", "group", "subtree_rate")), group, subtree_rate)

accel.nullscale=unique(accel.orig[,c("name", "null_scale")])

accel.new.wide<-spread(subset(accel.new, model=="neut_ver1", select=c("name", "group", "subtree_rate")), group, subtree_rate)

accel.rates<-merge(accel.new.wide, accel.nullscale, by="name", all.x=T, all.y=F)
accel.rates<-merge(accel.rates, accel.orig.wide, all.x=T, all.y=F, by="name", suffixes=c(".new", ".orig"))

#accel.tests<-read.table("moa/raw_convergence_results.tsv.gz", header=T, sep="\t")

filter.final=merge(filter.best, cnee.info, all=F)
filter.final<-merge(filter.final, accel.rates, all.x=T, all.y=T, by="name")
filter.final$clade.ct[is.na(filter.final$clade.ct)]=-1
