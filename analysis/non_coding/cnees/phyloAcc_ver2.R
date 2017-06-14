#new CNEE analysis
library(tidyr)
library(dplyr)
library(ggplot2)
setwd("/Volumes/LaCie/Projects/Current/ratites/final/phyloAcc_May2017/")

#functions
trunc_to_one <- function(x) {
  ifelse(x > 1, 1, x)
}


#read in data and clean up
lik<-read.table("Combined_elem_lik_061.txt", header=T, stringsAsFactors = F) %>% tbl_df %>%
  mutate(bf1 = loglik_RES - loglik_NUll, bf2 = loglik_RES - loglik_all) %>%
  select(cnee = ID, bf1, bf2, loglik.null = loglik_NUll, loglik.ratite = loglik_RES, loglik.full = loglik_all)

zpost<-read.table("Combined_post_Z_061.txt", header=T, stringsAsFactors = F) %>% tbl_df %>%
  mutate(cnee = lik$cnee)

zmat<-read.table("Combined_max_Z_061.txt", header=T, stringsAsFactors = F) %>% tbl_df

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
cnee<-inner_join(postacc, lik, by=c("cnee" = "cnee")) %>% inner_join(., postaccmat, by=c("cnee" = "cnee"), suffix=c(".prob", ".mat")) %>% ungroup

#define ratite accel classes
accel <- cnee %>% 
  mutate(rar5 = (bf1 > 5), rar10 = (bf1 > 10)) %>%
  mutate(rs5 = (bf2 > 5)) %>%
  mutate(conv1 = (ratite_loss_cons.prob >= 2), conv2 = (ratite_loss_cons_min.prob >= 2), conv3 = (ratite_loss_cons.mat >= 2), conv4 = (ratite_loss_cons_min.mat >= 2)) %>% select(cnee,  contains("ratite_loss"), rar5:conv4) 

#define sets
accel <- accel %>%
  mutate(set1 = rar5 & rs5, set2 = rar10 & rs5, set3 = rar5 & rs5 & nonratite_loss.mat == 0 & ratite_loss.mat >= 1, set4 = rar10 & rs5 & nonratite_loss.mat == 0 & ratite_loss.mat >= 1)

#use:
#set1, 3 <- 5 bf cutoff, 3 requires at least 1 ratite with post prob > 0.95 and no non-ratites with post prob > 0.95
#set2, 4 <- 10 bf cutoff, 4 requires at least 1 ratite with post prob > 0.95 and no non-ratites with post prob > 0.95
#conv1, conv2 <- biogeo and dollo using post probs
#conv3, conv4 <- biogeo and dollo using binary

#16 total comps: set1 - set4 * conv2 - conv5

write.table(accel, file="final_phyloAcc_cand_list.tsv", sep="\t", col.names = TRUE, row.names = FALSE)


##FUNCTIONAL ENRICHMENT###
#read in cnee <-> gene files

annot4<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/annotation/cnee/cnees.galgal4.annotation", sep="\t", stringsAsFactors = FALSE) %>% tbl_df %>% separate(V2, into=c("geneid", "symbol"), extra="merge")

annot5<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/annotation/cnee/cnees.galgal5.annotation", stringsAsFactors = FALSE, sep="\t") %>% tbl_df %>% separate(V2, into=c("geneid", "symbol"), extra="merge")


#merge
accel.gg4 <- full_join(accel, annot4, by=c("cnee" = "V1"))

##MESSING AROUND, NEEDS TO BE CLEANED UP##
#basic story very consistent with previous work, however

library(clusterProfiler)
library(DOSE)
foreground <- accel.gg4 %>% filter(set4==TRUE, conv2==TRUE) %>% distinct(geneid)
background <- accel.gg4 %>% distinct(geneid)

test1<-enrichGO(foreground$geneid, organism = "chicken", qvalueCutoff = 0.01, ont="BP")
test2<-compareCluster(foreground$geneid,organism="chicken")
summary(test1)
summary(test2)
?enrichKEGG
enrichMap(test1)
plotGOgraph(test1, firstSigNodes = 10)


library(PANTHER.db)
pthOrganisms(PANTHER.db)<-"CHICKEN"
PANTHER.db
columns(PANTHER.db)
cols<-c("PATHWAY_ID", "ENTREZ", "PATHWAY_TERM")
res<-select(PANTHER.db, keys=background$geneid, cols, "ENTREZ")
t2g <- data.frame(term=res$PATHWAY_ID, gene=res$ENTREZ)
t2n <- data.frame(term=res$PATHWAY_ID, name=res$PATHWAY_TERM)

test3<-enricher(foreground$geneid, universe = background$geneid, TERM2GENE = t2g, TERM2NAME = t2n, qvalueCutoff = 1, pvalueCutoff = 1)
summary(test3)
