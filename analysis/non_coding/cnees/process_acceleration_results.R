#process annotation results and make a raw annotation output file

library(data.table)
library(plyr)
accel.classes<-c("Casuar", "Rhea", "Kiwi", "Tinamou", "Ratite", "strCam", "anoDid", "BasalPaleo", "BasalPaleo_noR")
neut.mods<-c("neut_ver1", "neut_ver2", "neut_ver3")
run.types<-c("origAlign", "rrOnlyCNEEs", "withMoa")
base.path<-"/Volumes/LaCie/Projects/Current/ratites/final/accelTests/"

accel.res<-list()
for (run.type in run.types) {
  for (neut in neut.mods) {
    for (group in accel.classes) {
      print(paste0(c("Working on", run.type, neut, group), collapse=" "))
      
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
            
      file.type1<-paste0(base.path, run.type, "/", neut, "/", "all_", filegroup1, ".out", collapse="")
      file.type2<-paste0(base.path, run.type, "/", filegroup2, ".out_", neut, ".results", collapse="")
#debugging print statements
      print(file.type1)
      print(file.type2)
      if (file.exists(file.type1)) {
        infile<-file.type1
      } else {
        if (file.exists(file.type2)) {
          infile<-file.type2
        } else {
          next
        }
      }
      
      if (group=="BasalPaleo_noR" & neut != "neut_ver2") {
        next
      }
      
      #rr only run with neut2
      if (run.type=="rrOnlyCNEEs" & neut != "neut_ver2") {
        next
      }
       
      print("Got data")
      accel.res[[run.type]][[neut]][[group]]<-fread(infile, header=T, sep="\t", stringsAsFactors=F)
      accel.res[[run.type]][[neut]][[group]]$runtype = run.type
      accel.res[[run.type]][[neut]][[group]]$group = group
      accel.res[[run.type]][[neut]][[group]]$model = neut    
    }
  }
}


nm1.orig<-rbindlist(accel.res$orig$neut_ver1)
nm2.orig<-rbindlist(accel.res$orig$neut_ver2)
nm3.orig<-rbindlist(accel.res$orig$neut_ver3)
names(nm1.orig)[1]="chr"
names(nm2.orig)[1]="chr"
names(nm3.orig)[1]="chr"
all.rr_only<-rbindlist(accel.res$rrOnlyCNEEs$neut_ver2)
names(all.rr_only)[1]="chr"
nm1.moa<-rbindlist(accel.res$withMoa$neut_ver1)
nm2.moa<-rbindlist(accel.res$withMoa$neut_ver2)
nm3.moa<-rbindlist(accel.res$withMoa$neut_ver3)

#remove chr, start, end from all.orig and all.rr_only
nm1.orig[,chr:=NULL]
nm1.orig[,start:=NULL]
nm1.orig[,end:=NULL]
nm2.orig[,chr:=NULL]
nm2.orig[,start:=NULL]
nm2.orig[,end:=NULL]
nm3.orig[,chr:=NULL]
nm3.orig[,start:=NULL]
nm3.orig[,end:=NULL]
all.rr_only[,chr:=NULL]
all.rr_only[,start:=NULL]
all.rr_only[,end:=NULL]

#merge all into one data.table
all.accel<-rbind(nm1.orig, nm2.orig, nm3.orig, all.rr_only, nm1.moa, nm2.moa, nm3.moa)

save(all.accel, file=paste0(base.path, "all.accel.processed.Rdata", collapse=""))
