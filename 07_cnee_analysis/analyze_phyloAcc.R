#new CNEE analysis
library(tidyr)
library(dplyr)
library(ggplot2)
setwd("~/Projects/birds/ratite_compgen/data/phyloAcc_May2017/")

#functions
trunc_to_one <- function(x) {
  ifelse(x > 1, 1, x)
}


#read in data and clean up
lik<-read.table("Combined_elem_lik2.txt", header=T, stringsAsFactors = F) %>% tbl_df %>%
  mutate(bf1 = loglik_RES - loglik_NUll, bf2 = loglik_RES - loglik_all) %>%
  select(cnee = ID, bf1, bf2, loglik.null = loglik_NUll, loglik.ratite = loglik_RES, loglik.full = loglik_all)

#compare to old run
#lik_orig<-read.table("Combined_elem_lik_061.txt", header=T, stringsAsFactors = F) %>% tbl_df %>%
#  mutate(bf1 = loglik_RES - loglik_NUll, bf2 = loglik_RES - loglik_all) %>%
#  dplyr::select(cnee = ID, bf1, bf2, loglik.null = loglik_NUll, loglik.ratite = loglik_RES, loglik.full = #loglik_all)
#lik_comp <- inner_join(lik, lik_orig, by=c("cnee" = "cnee"), suffix = c(".new", ".orig"))
#ggplot(lik_comp, aes(x=bf1.new, y=bf1.orig)) + stat_binhex(bins=100)
#ggplot(lik_comp, aes(x=bf2.new, y=bf2.orig)) + stat_binhex(bins=100)
#lik_comp %>% filter(bf2.new <= 0, bf2.orig > 10)

zpost<-read.table("Combined_post_Z_06_2.txt", header=T, stringsAsFactors = F) %>% tbl_df %>%
  mutate(cnee = lik$cnee)

zmat<-read.table("Combined_max_Z_06_2.txt", header=T, stringsAsFactors = F) %>% tbl_df

#posterior prob acceleration
postacc<-zpost %>% select(contains("_3"))
colnames(postacc)<-colnames(zmat)
postacc$cnee = zpost$cnee

#convert to 1/0  matrix
postaccmat <- postacc[1:85]
postaccmat[postaccmat >= 0.95] = 1
postaccmat[postaccmat < 0.95] = 0
postaccmat$cnee = postacc$cnee

#posterior prob conservation
postcons<-zpost %>% select(contains("_2"))
colnames(postcons)<-colnames(zmat)
postcons$cnee = zpost$cnee

#define tips
vocal.tips<-c("taeGut", "ficAlb", "pseHum", "corBra", "melUnd", "calAnn")
ratite.tips<-c("aptRow", "aptOwe", "aptHaa", "casCas", "droNov", "rheAme", "rhePen", "anoDid", "strCam")
nonratite.tips <- colnames(zmat)[!colnames(zmat) %in% ratite.tips][1:27]
nonvl.tips <- colnames(zmat)[!colnames(zmat) %in% vocal.tips][1:30]

#compute expected number of losses for each ratite clade:
postacc <- postacc %>% mutate(cas_loss = (casCas - casCas.droNov) + (droNov - casCas.droNov) + (casCas.droNov - aptHaa.casCas)) %>%
  mutate(ost_loss = strCam - aptHaa.strCam) %>%
  mutate(rhea_loss = (rheAme - rheAme.rhePen) + (rhePen - rheAme.rhePen) + (rheAme.rhePen - aptHaa.rheAme)) %>%
  mutate(moa_loss = anoDid - cryCin.anoDid) %>%
  mutate(kiwi_loss = (aptHaa - aptHaa.aptRow) + (aptOwe - aptHaa.aptOwe) + (aptRow - aptHaa.aptRow) + (aptHaa.aptRow - aptHaa.aptOwe) + (aptHaa.aptOwe - aptHaa.casCas)) %>% 
  mutate(tin_loss = (tinGut - cryCin.tinGut) + (cryCin - cryCin.tinGut) + (eudEle - eudEle.notPer) + (notPer - eudEle.notPer) + (eudEle.notPer - cryCin.eudEle) + (cryCin.tinGut - cryCin.eudEle)) %>%
  mutate(internal_loss = (cryCin.eudEle - cryCin.anoDid) + (cryCin.anoDid - aptHaa.cryCin) + (aptHaa.rheAme - aptHaa.cryCin) + (aptHaa.casCas - aptHaa.rheAme) + (aptHaa.cryCin - aptHaa.strCam)) %>%
  mutate(ratite_loss = cas_loss + ost_loss + rhea_loss + moa_loss + kiwi_loss) %>%
  mutate(ratite_loss_cons = trunc_to_one(cas_loss) + trunc_to_one(rhea_loss) + trunc_to_one(kiwi_loss) + moa_loss + ost_loss) %>%
  mutate(ratite_loss_cons_min = trunc_to_one(cas_loss + rhea_loss + kiwi_loss) + moa_loss + ost_loss) %>%
    mutate(nonratite_loss = tin_loss + internal_loss)

postaccmat <- postaccmat %>% mutate(cas_loss = (casCas - casCas.droNov) + (droNov - casCas.droNov) + (casCas.droNov - aptHaa.casCas)) %>%
  mutate(ost_loss = strCam - aptHaa.strCam) %>%
  mutate(rhea_loss = (rheAme - rheAme.rhePen) + (rhePen - rheAme.rhePen) + (rheAme.rhePen - aptHaa.rheAme)) %>%
  mutate(moa_loss = anoDid - cryCin.anoDid) %>%
  mutate(kiwi_loss = (aptHaa - aptHaa.aptRow) + (aptOwe - aptHaa.aptOwe) + (aptRow - aptHaa.aptRow) + (aptHaa.aptRow - aptHaa.aptOwe) + (aptHaa.aptOwe - aptHaa.casCas)) %>% 
  mutate(tin_loss = (tinGut - cryCin.tinGut) + (cryCin - cryCin.tinGut) + (eudEle - eudEle.notPer) + (notPer - eudEle.notPer) + (eudEle.notPer - cryCin.eudEle) + (cryCin.tinGut - cryCin.eudEle)) %>%
  mutate(internal_loss = (cryCin.eudEle - cryCin.anoDid) + (cryCin.anoDid - aptHaa.cryCin) + (aptHaa.rheAme - aptHaa.cryCin) + (aptHaa.casCas - aptHaa.rheAme) + (aptHaa.cryCin - aptHaa.strCam)) %>%
  mutate(ratite_loss = cas_loss + ost_loss + rhea_loss + moa_loss + kiwi_loss) %>%
  mutate(ratite_loss_cons = trunc_to_one(cas_loss) + trunc_to_one(rhea_loss) + trunc_to_one(kiwi_loss) + moa_loss + ost_loss) %>%
  mutate(ratite_loss_cons_min = trunc_to_one(cas_loss + rhea_loss + kiwi_loss) + moa_loss + ost_loss) %>%
  mutate(nonratite_loss = tin_loss + internal_loss)

#compute expected neognath losses:
postacc <- postacc %>% mutate(neo_loss = (taeGut + ficAlb + pseHum + corBra + melUnd + falPer + picPub + lepDis + halLeu + aptFor + pygAde + fulGla + nipNip + balReg + chaVoc + calAnn + chaPel + cucCan +colLiv + mesUni + galGal + melGal + anaPla) - (taeGut.ficAlb + taeGut.pseHum + taeGut.corBra + taeGut.melUnd + taeGut.falPer + picPub.lepDis + picPub.halLeu + taeGut.picPub + aptFor.pygAde + aptFor.fulGla + aptFor.nipNip + taeGut.aptFor + balReg.chaVoc + taeGut.balReg + calAnn.chaPel + calAnn.cucCan + taeGut.calAnn + colLiv.mesUni + taeGut.colLiv + galGal.melGal + galGal.anaPla) - (2 * taeGut.galGal))

postaccmat <- postaccmat %>% mutate(neo_loss = (taeGut + ficAlb + pseHum + corBra + melUnd + falPer + picPub + lepDis + halLeu + aptFor + pygAde + fulGla + nipNip + balReg + chaVoc + calAnn + chaPel + cucCan +colLiv + mesUni + galGal + melGal + anaPla) - (taeGut.ficAlb + taeGut.pseHum + taeGut.corBra + taeGut.melUnd + taeGut.falPer + picPub.lepDis + picPub.halLeu + taeGut.picPub + aptFor.pygAde + aptFor.fulGla + aptFor.nipNip + taeGut.aptFor + balReg.chaVoc + taeGut.balReg + calAnn.chaPel + calAnn.cucCan + taeGut.calAnn + colLiv.mesUni + taeGut.colLiv + galGal.melGal + galGal.anaPla) - (2 * taeGut.galGal))

#make analysis dataset
phyloacc<-inner_join(postacc, lik, by=c("cnee" = "cnee")) %>% inner_join(., postaccmat, by=c("cnee" = "cnee"), suffix=c(".prob", ".mat")) %>% ungroup