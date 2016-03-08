
#process acceleration results

#load info
accel.classes<-c("Casuar", "Rhea", "Kiwi", "tinamou", "allRatite", "aptOwe", "aptHaa", "aptRow", "casCas", "droNov", "rheAme", "rhePen", "strCam", "anoDid", "basalPaleo")
neut.mods<-c("neut_ver1", "neut_ver2", "neut_ver3")
accel.res<-list()
for (neut in neut.mods) {
  for (group in accel.classes) {
    accel.res[[neut]][[group]]<-read.table(paste0(group, ".out_", neut, ".results", sep=""), header=T, sep="\t", comment.char="")
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

#moa.final has the moa acceleration results for each CNEE in long form
#now make a wide version that has one row per CNEE

library(tidyr)
moa.sub<-subset(moa.final, select=c("name", "qval", "group", "model"))
moa.wide<-spread(moa.sub, group, qval)

#test different neutral models
accel.sensitivity<-data.frame(mod1=character(0), mod2=character(0), group=character(0), ff=numeric(0), tt=numeric(0), tf=numeric(0), total=numeric(0))
for (neut in neut.mods) {
	for (neut2 in neut.mods) {
		for (group in accel.classes) {
			temp.table<-table(moa.wide[moa.wide$model==neut,group]<0.01, moa.wide[moa.wide$model==neut2,group]<0.01)
			ff<-temp.table[1,1]
			tf<-temp.table[1,2]+temp.table[2,1]
			tt<-temp.table[2,2]
			total<-ff+tf+tt
			temp.df<-data.frame(mod1=neut, mod2=neut2, group=group, ff=ff, tt=tt, tf=tf, total=total)
			accel.sensitivity<-rbind(accel.sensitivity, temp.df)
		}
	}
}

accel.sensitivity$comp=apply(accel.sensitivity[,c(1,2)], 1, paste0, collapse=":")
accel.sensitivity$prop.diff=accel.sensitivity$tf/(accel.sensitivity$tt+accel.sensitivity$tf)

#add code to write out moa.wide data frame

