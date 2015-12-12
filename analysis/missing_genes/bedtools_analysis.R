setwd("~/Desktop/maker/out")
path=("~/Desktop/maker/out")

file.names <- dir(path, pattern =".txt")

#read in all files that end in .txt from maker/out and place them into one dataframe
best.out<-""
for(i in 1:length(file.names)){
  file <- read.delim(file.names[i])
  best.out <- Reduce(function(x, y) merge(x, y, all=TRUE), list(best.out,file))
}

drops <- c("x") #create a list of the columns that will be uninformative soon and we'd like to drop (just the extra x)
nbest.out <- best.out[,!(names(best.out) %in% drops)]

#_Lift values == number of exons lifted over from chicken in a given species
#_Inter values == number of lifted over exons that intersect (even for just one base) with an exon from a gene model in that species

#some sanity checks: chicken should have equal numbers in its _Lift and _Inter values
all(nbest.out$galGal_bed12_Lift==nbest.out$galGal_bed12_Inter) #TRUE
all_animals <- c("chrPic_bed12_Inter","chrPic_bed12_Lift","allMis_bed12_Inter","allMis_bed12_Lift","anoCar_bed12_Inter","anoCar_bed12_Lift","aptHaa_bed12_Inter","aptOwe_bed12_Inter","aptRow_bed12_Inter","casCas_bed12_Inter","droNov_bed12_Inter","rheAme_bed12_Inter","rhePen_bed12_Inter","strCam_bed12_Inter","aptHaa_bed12_Lift","aptOwe_bed12_Lift","aptRow_bed12_Lift","casCas_bed12_Lift","droNov_bed12_Lift","rheAme_bed12_Lift","rhePen_bed12_Lift","strCam_bed12_Lift","galGal_bed12_Inter","anaPla_bed12_Inter","aptFor_bed12_Inter","balReg_bed12_Inter","calAnn_bed12_Inter","chaPel_bed12_Inter","chaVoc_bed12_Inter","colLiv_bed12_Inter","corBra_bed12_Inter","cryCin_bed12_Inter","cucCan_bed12_Inter","eudEle_bed12_Inter","falPer_bed12_Inter","ficAlb_bed12_Inter","fulGla_bed12_Inter","halLeu_bed12_Inter","lepDis_bed12_Inter","melGal_bed12_Inter","melUnd_bed12_Inter","mesUni_bed12_Inter","nipNip_bed12_Inter","notPer_bed12_Inter","picPub_bed12_Inter","pseHum_bed12_Inter","pygAde_bed12_Inter","taeGut_bed12_Inter","tinGut_bed12_Inter","galGal_bed12_Lift","anaPla_bed12_Lift","aptFor_bed12_Lift","balReg_bed12_Lift","calAnn_bed12_Lift","chaPel_bed12_Lift","chaVoc_bed12_Lift","colLiv_bed12_Lift","corBra_bed12_Lift","cryCin_bed12_Lift","cucCan_bed12_Lift","eudEle_bed12_Lift","falPer_bed12_Lift","ficAlb_bed12_Lift","fulGla_bed12_Lift","halLeu_bed12_Lift","lepDis_bed12_Lift","melGal_bed12_Lift","melUnd_bed12_Lift","mesUni_bed12_Lift","nipNip_bed12_Lift","notPer_bed12_Lift","picPub_bed12_Lift","pseHum_bed12_Lift","pygAde_bed12_Lift","taeGut_bed12_Lift","tinGut_bed12_Lift")
nbest.out$max_test <- apply(nbest.out[,all_animals], 1, function(x) max(x, na.rm=T))
all(nbest.out$galGal_bed12_Lift==nbest.out$max_test)

#drop sanity check "max_test" column
drops <- c("max_test")
nbest.out <- nbest.out[,!(names(nbest.out) %in% drops)]

#tim and phil first pass code with names changed after re-generation of data (straight from BED files instead of converting from PSL)
ratite_inter <- c("aptHaa_bed12_Inter","aptOwe_bed12_Inter","aptRow_bed12_Inter","casCas_bed12_Inter","droNov_bed12_Inter","rheAme_bed12_Inter","rhePen_bed12_Inter","strCam_bed12_Inter")
ratite_lift <- c("aptHaa_bed12_Lift","aptOwe_bed12_Lift","aptRow_bed12_Lift","casCas_bed12_Lift","droNov_bed12_Lift","rheAme_bed12_Lift","rhePen_bed12_Lift","strCam_bed12_Lift")
flying_inter <-c("galGal_bed12_Inter","anaPla_bed12_Inter","aptFor_bed12_Inter","balReg_bed12_Inter","calAnn_bed12_Inter","chaPel_bed12_Inter","chaVoc_bed12_Inter","colLiv_bed12_Inter","corBra_bed12_Inter","cryCin_bed12_Inter","cucCan_bed12_Inter","eudEle_bed12_Inter","falPer_bed12_Inter","ficAlb_bed12_Inter","fulGla_bed12_Inter","halLeu_bed12_Inter","lepDis_bed12_Inter","melGal_bed12_Inter","melUnd_bed12_Inter","mesUni_bed12_Inter","nipNip_bed12_Inter","notPer_bed12_Inter","picPub_bed12_Inter","pseHum_bed12_Inter","pygAde_bed12_Inter","taeGut_bed12_Inter","tinGut_bed12_Inter")
flying_lift <-c("galGal_bed12_Lift","anaPla_bed12_Lift","aptFor_bed12_Lift","balReg_bed12_Lift","calAnn_bed12_Lift","chaPel_bed12_Lift","chaVoc_bed12_Lift","colLiv_bed12_Lift","corBra_bed12_Lift","cryCin_bed12_Lift","cucCan_bed12_Lift","eudEle_bed12_Lift","falPer_bed12_Lift","ficAlb_bed12_Lift","fulGla_bed12_Lift","halLeu_bed12_Lift","lepDis_bed12_Lift","melGal_bed12_Lift","melUnd_bed12_Lift","mesUni_bed12_Lift","nipNip_bed12_Lift","notPer_bed12_Lift","picPub_bed12_Lift","pseHum_bed12_Lift","pygAde_bed12_Lift","taeGut_bed12_Lift","tinGut_bed12_Lift")
tinamou_inter <- c("cryCin_bed12_Inter","eudEle_bed12_Inter","notPer_bed12_Inter","tinGut_bed12_Inter")
tinamou_lift <- c("cryCin_bed12_Lift","eudEle_bed12_Lift","notPer_bed12_Lift","tinGut_bed12_Lift")


#don't worry about errors.  they crop up for ratites and tinamous if all species have NA at specific gene ID.
nbest.out$r_lowLift <- apply(nbest.out[,ratite_lift], 1, function(x) min(x, na.rm=T))
nbest.out$f_lowLift <- apply(nbest.out[,flying_lift], 1, function(x) min(x, na.rm=T))
nbest.out$t_lowLift <- apply(nbest.out[,tinamou_lift], 1, function(x) min(x, na.rm=T))

nbest.out$r_lowInt <- apply(nbest.out[,ratite_inter], 1, function(x) min(x, na.rm=T))
nbest.out$f_lowInt <- apply(nbest.out[,flying_inter], 1, function(x) min(x, na.rm=T))
nbest.out$t_lowInt <- apply(nbest.out[,tinamou_inter], 1, function(x) min(x, na.rm=T))

nbest.out$r_hiInt <- apply(nbest.out[,ratite_inter], 1, function(x) max(x, na.rm=T))

nbest.out$ratite_lost = nbest.out$r_hiInt == 0 & nbest.out$r_lowLift > 0 & nbest.out$r_lowLift < Inf
 
ratite.lowfruit <- subset(nbest.out, nbest.out$ratite_lost & nbest.out$t_lowInt > 0 & nbest.out$t_lowInt < Inf)
ratite.lowfruit2 <- subset(nbest.out, nbest.out$ratite_lost & nbest.out$t_lowInt > 0 & nbest.out$t_lowInt < Inf)

write.csv(nbest.out, "liftover_gene_model_exon_intersections.csv", row.names = FALSE)