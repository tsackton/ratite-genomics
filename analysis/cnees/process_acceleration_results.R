#process annotation results and make a raw annotation output file

library(data.table)
library(plyr)
accel.classes<-c("Casuar", "Rhea", "Kiwi", "Tinamou", "Ratite", "aptOwe", "aptHaa", "aptRow", "casCas", "droNov", "rheAme", "rhePen", "strCam", "anoDid", "BasalPaleo", "BasalPaleo_noR")
neut.mods<-c("neut_ver1", "neut_ver2", "neut_ver3")
run.types<-c("orig", "rrOnlyCNEEs", "withMoa")
base.path<-"/Volumes/LaCie/Projects/Current/ratites/final/accelTests/"

accel.res<-list()
for (run.type in run.types) {
  for (neut in neut.mods) {
    for (group in accel.classes) {
      #complicated code to deal with my annoyingly different file naming conventions
      filegroup1=group
      filegroup2=group
      if (group == "Ratite") {
        filegroup1 = tolower(filegroup1)
        filegroup2 = "allRatite"
      }
      if (group == "Tinamou") {
        filegroup1 = tolower(filegroup1)
      }
      if (group == "BasalPaleo") {
        filegroup1 = "basalPaleo"
        filegroup2 = "basalPaleo"
      }
      if (group == "BasalPaleo_noR") {
        filegroup1 = "norr_basalPaleo"
      }
      
      file.type1<-paste0(base.path, run.type, "/", "all_", filegroup1, ".out", collapse="")
      file.type2<-paste0(base.path, run.type, "/", filegroup2, ".out_", neut, ".results", collapse="")
      if (file.exists(file.type1)) {
        infile<-file.type1
      } else {
        if (file.exists(file.type2)) {
          infile<-file.type2
        } else {
          next
        }
      }
      
      #orig and rr only run with neut2
      if (((run.type=="orig")| (run.type=="rrOnlyCNEEs")) & neut != "neut_ver2") {
        next
      }
      
      print(paste0(c("Working on", run.type, neut, group), collapse=" "))
      
      accel.res[[run.type]][[neut]][[group]]<-fread(infile, header=T, sep="\t", stringsAsFactors=F)
      accel.res[[run.type]][[neut]][[group]]$runtype = run.type
      accel.res[[run.type]][[neut]][[group]]$group = group
      accel.res[[run.type]][[neut]][[group]]$model = neut    
    }
  }
}


all.orig<-rbindlist(accel.res$orig$neut_ver2)
names(all.orig)[1]="chr"
all.rr_only<-rbindlist(accel.res$rrOnlyCNEEs$neut_ver2)
names(all.rr_only)[1]="chr"
nm1.moa<-rbindlist(accel.res$withMoa$neut_ver1)
nm2.moa<-rbindlist(accel.res$withMoa$neut_ver2)
nm3.moa<-rbindlist(accel.res$withMoa$neut_ver3)

#remove chr, start, end from all.orig and all.rr_only
all.orig[,chr:=NULL]
all.orig[,start:=NULL]
all.orig[,end:=NULL]
all.rr_only[,chr:=NULL]
all.rr_only[,start:=NULL]
all.rr_only[,end:=NULL]

#merge all into one data.table
all.accel<-rbind(all.orig, all.rr_only, nm1.moa, nm2.moa, nm3.moa)

save(all.accel, file=paste0(base.path, "all.accel.processed.Rdata", collapse=""))
