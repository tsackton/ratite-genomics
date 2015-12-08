setwd("~/Desktop/maker/out")
path=("~/Desktop/maker/out")

file.names <- dir(path, pattern =".txt")

best.out<-""
for(i in 1:length(file.names)){
  file <- read.delim(file.names[i])
  best.out <- Reduce(function(x, y) merge(x, y, all=TRUE), list(best.out,file))
}

drops <- c("x") #create a list of the columns that will be uninformative soon and we'd like to drop (just the extra x)
nbest.out <- best.out[,!(names(best.out) %in% drops)]

#some sanity checks
all(nbest.out$galGal_bed12_Lift==nbest.out$galGal_bed12_Inter)
all_animals <- c("chrPic_bed12_Inter","chrPic_bed12_Lift","allMis_bed12_Inter","allMis_bed12_Lift","anoCar_bed12_Inter","anoCar_bed12_Lift","aptHaa_bed12_Inter","aptOwe_bed12_Inter","aptRow_bed12_Inter","casCas_bed12_Inter","droNov_bed12_Inter","rheAme_bed12_Inter","rhePen_bed12_Inter","strCam_bed12_Inter","aptHaa_bed12_Lift","aptOwe_bed12_Lift","aptRow_bed12_Lift","casCas_bed12_Lift","droNov_bed12_Lift","rheAme_bed12_Lift","rhePen_bed12_Lift","strCam_bed12_Lift","galGal_bed12_Inter","anaPla_bed12_Inter","aptFor_bed12_Inter","balReg_bed12_Inter","calAnn_bed12_Inter","chaPel_bed12_Inter","chaVoc_bed12_Inter","colLiv_bed12_Inter","corBra_bed12_Inter","cryCin_bed12_Inter","cucCan_bed12_Inter","eudEle_bed12_Inter","falPer_bed12_Inter","ficAlb_bed12_Inter","fulGla_bed12_Inter","halLeu_bed12_Inter","lepDis_bed12_Inter","melGal_bed12_Inter","melUnd_bed12_Inter","mesUni_bed12_Inter","nipNip_bed12_Inter","notPer_bed12_Inter","picPub_bed12_Inter","pseHum_bed12_Inter","pygAde_bed12_Inter","taeGut_bed12_Inter","tinGut_bed12_Inter","galGal_bed12_Lift","anaPla_bed12_Lift","aptFor_bed12_Lift","balReg_bed12_Lift","calAnn_bed12_Lift","chaPel_bed12_Lift","chaVoc_bed12_Lift","colLiv_bed12_Lift","corBra_bed12_Lift","cryCin_bed12_Lift","cucCan_bed12_Lift","eudEle_bed12_Lift","falPer_bed12_Lift","ficAlb_bed12_Lift","fulGla_bed12_Lift","halLeu_bed12_Lift","lepDis_bed12_Lift","melGal_bed12_Lift","melUnd_bed12_Lift","mesUni_bed12_Lift","nipNip_bed12_Lift","notPer_bed12_Lift","picPub_bed12_Lift","pseHum_bed12_Lift","pygAde_bed12_Lift","taeGut_bed12_Lift","tinGut_bed12_Lift")
nbest.out$max_test <- apply(nbest.out[,all_animals], 1, function(x) max(x, na.rm=T))
all(nbest.out$galGal_bed12_Lift==nbest.out$max_test)


#tim and phil first pass
ratite_inter <- c("aptHaa_bed12_Inter","aptOwe_bed12_Inter","aptRow_bed12_Inter","casCas_bed12_Inter","droNov_bed12_Inter","rheAme_bed12_Inter","rhePen_bed12_Inter","strCam_bed12_Inter")
ratite_lift <- c("aptHaa_bed12_Lift","aptOwe_bed12_Lift","aptRow_bed12_Lift","casCas_bed12_Lift","droNov_bed12_Lift","rheAme_bed12_Lift","rhePen_bed12_Lift","strCam_bed12_Lift")
flying_inter <-c("galGal_bed12_Inter","anaPla_bed12_Inter","aptFor_bed12_Inter","balReg_bed12_Inter","calAnn_bed12_Inter","chaPel_bed12_Inter","chaVoc_bed12_Inter","colLiv_bed12_Inter","corBra_bed12_Inter","cryCin_bed12_Inter","cucCan_bed12_Inter","eudEle_bed12_Inter","falPer_bed12_Inter","ficAlb_bed12_Inter","fulGla_bed12_Inter","halLeu_bed12_Inter","lepDis_bed12_Inter","melGal_bed12_Inter","melUnd_bed12_Inter","mesUni_bed12_Inter","nipNip_bed12_Inter","notPer_bed12_Inter","picPub_bed12_Inter","pseHum_bed12_Inter","pygAde_bed12_Inter","taeGut_bed12_Inter","tinGut_bed12_Inter")
flying_lift <-c("galGal_bed12_Lift","anaPla_bed12_Lift","aptFor_bed12_Lift","balReg_bed12_Lift","calAnn_bed12_Lift","chaPel_bed12_Lift","chaVoc_bed12_Lift","colLiv_bed12_Lift","corBra_bed12_Lift","cryCin_bed12_Lift","cucCan_bed12_Lift","eudEle_bed12_Lift","falPer_bed12_Lift","ficAlb_bed12_Lift","fulGla_bed12_Lift","halLeu_bed12_Lift","lepDis_bed12_Lift","melGal_bed12_Lift","melUnd_bed12_Lift","mesUni_bed12_Lift","nipNip_bed12_Lift","notPer_bed12_Lift","picPub_bed12_Lift","pseHum_bed12_Lift","pygAde_bed12_Lift","taeGut_bed12_Lift","tinGut_bed12_Lift")
tinamou_inter <- c("cryCin_bed12_Inter","eudEle_bed12_Inter","notPer_bed12_Inter","tinGut_bed12_Inter")
tinamou_lift <- c("cryCin_bed12_Lift","eudEle_bed12_Lift","notPer_bed12_Lift","tinGut_bed12_Lift")

#MAKE NOTE 
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






#old code below:

# nbest.out[is.na(nbest.out)] <- 0 #change NA values to 0
# 
# 
# colnames(nbest.out) #used to list the column names in order to create the indexes below
# #below are the three indexes that we can use to filter through the data.  rindex is ratites, ntindex is neognaths and tinamous, rest are herps
# ratite_inter <- c("aptHaa_Inter","aptOwe_Inter","aptRow_Inter","casCas_Inter","droNov_Inter","rheAme_Inter","rhePen_Inter")
# ratite_lift <- c("aptHaa_Lift","aptOwe_Lift","aptRow_Lift","casCas_Lift","droNov_Lift","rheAme_Lift","rhePen_Lift")
# ratite_both <- c("aptHaa_Inter","aptOwe_Inter","aptRow_Inter","casCas_Inter","droNov_Inter","rheAme_Inter","rhePen_Inter","aptHaa_Lift","aptOwe_Lift","aptRow_Lift","casCas_Lift","droNov_Lift","rheAme_Lift","rhePen_Lift")
# 
# flying_inter <-c("galGal_Inter","anaPla_Inter","aptFor_Inter","balReg_Inter","calAnn_Inter","chaPel_Inter","chaVoc_Inter","colLiv_Inter","corBra_Inter","cryCin_Inter","cucCan_Inter","eudEle_Inter","falPer_Inter","ficAlb_Inter","fulGla_Inter","halLeu_Inter","lepDis_Inter","melGal_Inter","melUnd_Inter","mesUni_Inter","nipNip_Inter","notPer_Inter","picPub_Inter","pseHum_Inter","pygAde_Inter","taeGut_Inter","tinGut_Inter")
# flying_lift <-c("galGal_Lift","anaPla_Lift","aptFor_Lift","balReg_Lift","calAnn_Lift","chaPel_Lift","chaVoc_Lift","colLiv_Lift","corBra_Lift","cryCin_Lift","cucCan_Lift","eudEle_Lift","falPer_Lift","ficAlb_Lift","fulGla_Lift","halLeu_Lift","lepDis_Lift","melGal_Lift","melUnd_Lift","mesUni_Lift","nipNip_Lift","notPer_Lift","picPub_Lift","pseHum_Lift","pygAde_Lift","taeGut_Lift","tinGut_Lift")
# flying_both <- c("galGal_Inter","anaPla_Inter","aptFor_Inter","balReg_Inter","calAnn_Inter","chaPel_Inter","chaVoc_Inter","colLiv_Inter","corBra_Inter","cryCin_Inter","cucCan_Inter","eudEle_Inter","falPer_Inter","ficAlb_Inter","fulGla_Inter","halLeu_Inter","lepDis_Inter","melGal_Inter","melUnd_Inter","mesUni_Inter","nipNip_Inter","notPer_Inter","picPub_Inter","pseHum_Inter","pygAde_Inter","taeGut_Inter","tinGut_Inter","galGal_Lift","anaPla_Lift","aptFor_Lift","balReg_Lift","calAnn_Lift","chaPel_Lift","chaVoc_Lift","colLiv_Lift","corBra_Lift","cryCin_Lift","cucCan_Lift","eudEle_Lift","falPer_Lift","ficAlb_Lift","fulGla_Lift","halLeu_Lift","lepDis_Lift","melGal_Lift","melUnd_Lift","mesUni_Lift","nipNip_Lift","notPer_Lift","picPub_Lift","pseHum_Lift","pygAde_Lift","taeGut_Lift","tinGut_Lift")
# 
# tinamou_inter <- c("cryCin_Inter","eudEle_Inter","notPer_Inter","tinGut_Inter")
# tinamou_lift <- c("cryCin_Lift","eudEle_Lift","notPer_Lift","tinGut_Lift")
# tinamou_both <- c("cryCin_Inter","eudEle_Inter","notPer_Inter","tinGut_Inter","cryCin_Lift","eudEle_Lift","notPer_Lift","tinGut_Lift")
# 
# paleo_inter <- c(ratite_inter,tinamou_inter)
# paleo_lift <- c(ratite_lift,tinamou_lift)
# paleo_both <- c(ratite_both,tinamou_both)
# 
# herps_inter <- c("anoCar_Inter","allMis_Inter","chrPic_Inter")
# herps_lift <- c("anoCar_Lift","allMis_Lift","chrPic_Lift")
# herps_both <- c("anoCar_Inter","allMis_Inter","chrPic_Inter","anoCar_Lift","allMis_Lift","chrPic_Lift")
# 
# #this was the old way, but was not versatile
# #not updated ratite.nlowfruit <- subset(best.out, best.out$rheAme > 1 & best.out$rhePen > 1 & best.out$casCas > 1 & best.out$droNov > 1 & best.out$strCam > 1 & best.out$aptHaa > 1 & best.out$aptOwe > 1 & best.out$aptRow > 1 & best.out$eudEle <= 1 & best.out$cryCin <= 1 & best.out$notPer <= 1 & best.out$melUnd <= 1 & best.out$pseHum <= 1 & best.out$colLiv <= 1 & best.out$falPer <= 1 & best.out$anaPla <= 1 & best.out$fulGla <= 1 & best.out$lepDis <= 1 & best.out$corBra <= 1 & best.out$mesUni <= 1 & best.out$picPub <= 1 & best.out$calAnn <= 1 & best.out$pygAde <= 1 & best.out$aptFor <= 1 & best.out$chaVoc <= 1 & best.out$nipNip <= 1 & best.out$cucCan <= 1 & best.out$balReg <= 1 & best.out$halLeu <= 1 & best.out$chaPel & best.out$tinGut <= 1)
# 
# #want the subset of GeneIDs where ratites have same liftover number as chicken
# nbest.out$r_lowLift <- apply(nbest.out[,ratite_lift], 1, function(x) min(x))
# nbest.out$f_lowLift <- apply(nbest.out[,flying_lift], 1, function(x) min(x))
# 
# nbest.out$r_lowInt <- apply(nbest.out[,ratite_inter], 1, function(x) min(x))
# nbest.out$f_lowInt <- apply(nbest.out[,flying_inter], 1, function(x) min(x))
# 
# nbest.out$r_hiInt <- apply(nbest.out[,ratite_inter], 1, function(x) max(x))
# nbest.out$r_IntDec <- nbest.out$r_hiInt/nbest.out$galGal_Inter 
# 
# nbest.out$t_hiInt <- apply(nbest.out[,tinamou_inter], 1, function(x) max(x))
# nbest.out$t_lowInt <- apply(nbest.out[,tinamou_inter], 1, function(x) min(x))
# 
# nbest.out$f_hiInt <- apply(nbest.out[,flying_inter], 1, function(x) max(x))
# 
# nbest.out$sumIntFly <- apply(nbest.out[,flying_inter], 1, function(x) sum(x)/length(x))
# nbest.out$propFly <- nbest.out$sumIntFly/nbest.out$galGal_Inter
# 
# rat.sameLO <- subset(nbest.out, nbest.out$r_lowLift == best.out$galGal_Lift) 
# bird.sameLO <- subset(nbest.out, nbest.out$f_lowLift > 0  & nbest.out$r_lowLift > 0) #both liftover values are same as chicken across the birds
# 
# ratite.lowfruit <- subset(bird.sameLO,bird.sameLO$r_hiInt < bird.sameLO$galGal_Inter & bird.sameLO$f_lowInt == bird.sameLO$galGal_Inter)
# ratite.lowfruit2 <- subset(bird.sameLO, bird.sameLO$r_hiInt == 0 & bird.sameLO$t_lowInt == bird.sameLO$galGal_Inter & bird.sameLO$propFly > 0.5)
# 
# tim.firsttest <- subset(bird.sameLO, bird.sameLO$)
# 
# 
# 
# 
# drops <- c("f_IntDec") #create a list of the columns that will be uninformative soon and we'd like to drop (just the extra x)
# nbest.out <- nbest.out[,!(names(nbest.out) %in% drops)]
# 
# nbest.out$gGequalLO <- apply(nbest.out[,ratite_lift], 1, function(x) sum(x > 1)/length(x)) #create a new column called badRatites that contains the proportion of ratites with more than 1 stop for each gene
# nbest.out$goodFliers <- apply(nbest.out[,ntindexn], 1, function(x) sum(x <= 1)/length(x)) #create a new column called goodFliers that contains the proportion of ntindex birds with 1 or 0 stops for each gene
# nbest.out$present_nt <- apply(nbest.out[,ntindexn], 1, function(x) sum(x > 0)/length(x))
# nbest.out$goodTinamous <- apply(nbest.out[,tindexn], 1, function(x) sum(x <= 1)/length(x))
# nbest.out$present_palaeo <- apply(nbest.out[,pindexn], 1, function(x) sum(x > 0)/length(x))
# nbest.out$sumRatites <- apply(nbest.out[,rindexn], 1, function(x) sum(x > 5))
# nbest.out$absent_herps <- apply(nbest.out[,restn], 1, function(x) sum(x < 0)/length(x))
# 
# 
# ratite.nlowfruit <- subset(nbest.out, nbest.out$badRatites>0.8 & nbest.out$goodFliers>0.99) #example of how to filter and subset data using these - much cleaner and more versatile than above subsetting
# ratite.nlowfruit2 <- subset(nbest.out, nbest.out$badRatites>=0.8 & nbest.out$goodFliers>0.9 & nbest.out$genePresent>0.5) 
# ratite.nlowfruit3 <- subset(nbest.out, nbest.out$badRatites>=0.8 & nbest.out$goodFliers>0.8 & nbest.out$genePresent==1) 
# paleo.nfirstlook <- subset(nbest.out, nbest.out$badRatites==1 & nbest.out$goodTinamous == 1 & nbest.out$present_palaeo==1) # bad in all ratites, good in all tinamous, present in all paleo
# ratite.nlowfruit4 <- subset(nbest.out, nbest.out$badRatites>=0.5 & nbest.out$goodTinamous == 1 & nbest.out$present_palaeo==1 & nbest.out$present_nt>0.8 & nbest.out$goodFliers>0.99) # bad in some ratites, good and present in everything else
# ratite.nlowfruit5 <- subset(nbest.out, nbest.out$badRatites>0 & nbest.out$badRatites< 0.5 & nbest.out$goodTinamous == 1 & nbest.out$present_palaeo==1 & nbest.out$present_nt>0.8 & nbest.out$goodFliers>0.99) # bad in some ratites, good and present in everything else
# ratite.nlowfruit6 <- subset(nbest.out, nbest.out$sumRatites==8 & nbest.out$present_palaeo==1 & nbest.out$goodTinamous >= 0.75 & nbest.out$goodFliers >= 0.8 & nbest.out$present_nt >= 0.8)
# ratite.nlowfruit7 <- subset(nbest.out, nbest.out$goodFliers==1 & nbest.out$sumRatites>3 & nbest.out$present_nt > 0.8 & nbest.out$present_palaeo >0.8)
# 
