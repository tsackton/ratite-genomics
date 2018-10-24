## CODE TO PARSE PHYLOACC OUTPUT ##
## UPDATED OCT 2018 FOR MANUSCRIPT REVISIONS ##

library(tidyverse)

#functions
tto <- function(x) {
  ifelse(x > 1, 1, x)
}

ctb <- function(x, cutoff=0.90, lower=FALSE) {
  if (lower) {
    ifelse(x <= cutoff, 1, 0)
  }
  else {
    ifelse(x >= cutoff, 1, 0)
  }
}

get_max <- function(x, cutoff = 0.90) {
  ifelse(max(x) > cutoff, 1, 0)
}

input_lik <- function(file, dataset) {
  read_delim(file=file, delim="\t") %>%
    mutate(dataset = dataset) %>%
    rename(cnee = ID) %>%
    select(-No.)
}

input_states <- function(file, dataset, datatype) {
  #hacky but hopefully works
  
  #get column names:
  col_names <- read_delim(file, delim="\t", n_max=1, col_names=FALSE) %>% slice(1:1) %>% unlist() %>% unname %>% c("cnee", .)
  
  #col types depends on whether this is posterior probabilty matrix (double) or maxium posterior state (integer)
  if (datatype == "postprob") {
    col_types <- paste0(c("c", rep("d", length(col_names)-1)), collapse="")
  }
  else if (datatype == "maxpost") {
    col_types <- paste0(c("c", rep("i", length(col_names)-1)), collapse="")
  }
  else {
    stop("need to specify posterior probabilities (datatype=postprob) or maxium a posteriori state (datatyp=maxpost) with datatype argument")
  }
  read_delim(file=file, delim="\t", skip=1, col_names = col_names, col_types = col_types) %>% 
    mutate(dataset = dataset)
}

merge_phyloAcc <- function(dataset) {
  #merge different runs
  dfpost <- bind_rows(posteriors[[dataset]], .id="version") %>%
    select(cnee, version, dataset, n_rate, c_rate, g_rate, l_rate, ends_with("_3")) %>%
    rename_at(vars(ends_with("_3")), gsub, pattern="_3", replacement="", fixed=TRUE)
  dflik <- bind_rows(likelihoods[[dataset]], .id="version")
  
  #merge posteriors and likelihoods
  inner_join(dflik, dfpost, by=c("cnee" = "cnee", "version" = "version", "dataset" = "dataset")) %>% 
  select(dataset, version, cnee, everything())
}

## END FUNCTIONS ##

likelihoods <- list()
posteriors <- list()

#read in files
alignments<-c("EXTEND_1012", "ORG_1012", "REDUCE_1012")
versions<-c(".", "gain", "gain_gap", "gap")

for (align in alignments) {
  for (ver in versions) {
    path = paste0("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/03_phyloAcc_res/", align, "/", ver, "/")
    ver = ifelse(ver == ".", "orig", ver)
    dataset = paste0(align, "-", ver)
    likelihoods[[align]][[ver]] = input_lik(file=paste0(path, "Combined_elem_lik.txt"), dataset)
    posteriors[[align]][[ver]] = input_states(file=paste0(path, "Combined_post_Z_M2.txt"), dataset, "postprob")
  }
}

#now merge each dataset to make one cnee set for each phyloAcc run, containing the CNEE id, basic parameters, post prob acceleration, max a posteriori estimated state

cnee_orig <- merge_phyloAcc("ORG_1012")
cnee_reduced <- merge_phyloAcc("REDUCE_1012")
cnee_ext <- merge_phyloAcc("EXTEND_1012")

#clean up 
rm(likelihoods)
rm(posteriors)

#now, need to compute some stats on these merged datasets. in particular, want expected number of losses for each flightless clade plus total flightless
#also need tinamou losses and volant neognath losses (penguin separate)
#for each loss calculation, do the following:
#assume dollo (e.g. e/k/c/r just one clade) or assume biogeographic (e.g. e/c, k, r, o, m all separate)
#use cutoff vs use raw posterior probs

#this is going to be a bit of a mess, need to name columns carefully

cnee_orig_losses <- cnee_orig %>%
  #this is the posterior estimated number of losses on each clade
  mutate(cd_pp_loss = (casCas - casCas.droNov) + (droNov - casCas.droNov) + (casCas.droNov - aptHaa.casCas)) %>%
  mutate(rh_pp_loss = (rheAme - rheAme.rhePen) + (rhePen - rheAme.rhePen) + (rheAme.rhePen - aptHaa.rheAme)) %>%
  mutate(os_pp_loss = strCam - aptHaa.strCam) %>%
  mutate(ki_pp_loss = (aptHaa - aptHaa.aptRow) + (aptOwe - aptHaa.aptOwe) + (aptRow - aptHaa.aptRow) + (aptHaa.aptRow - aptHaa.aptOwe) + (aptHaa.aptOwe - aptHaa.casCas)) %>%
  mutate(mo_pp_loss = anoDid - cryCin.anoDid) %>%
  mutate(ti_pp_loss = (tinGut - cryCin.tinGut) + (cryCin - cryCin.tinGut) + (eudEle - eudEle.notPer) + (notPer - eudEle.notPer) + (eudEle.notPer - cryCin.eudEle) + (cryCin.tinGut - cryCin.eudEle)) %>%
  mutate(it_pp_loss = (cryCin.eudEle - cryCin.anoDid) + (cryCin.anoDid - aptHaa.cryCin) + (aptHaa.rheAme - aptHaa.cryCin) + (aptHaa.casCas - aptHaa.rheAme) + (aptHaa.cryCin - aptHaa.strCam)) %>%
  mutate(neo_pp_loss = (taeGut + ficAlb + pseHum + corBra + melUnd + falPer + picPub + lepDis + halLeu + fulGla + nipNip + balReg + chaVoc + calAnn + chaPel + cucCan +colLiv + mesUni + galGal + melGal + anaPla) - (taeGut.ficAlb + taeGut.pseHum + taeGut.corBra + taeGut.melUnd + taeGut.falPer + picPub.lepDis + picPub.halLeu + taeGut.picPub  + aptFor.nipNip + taeGut.aptFor + balReg.chaVoc + taeGut.balReg + calAnn.chaPel + calAnn.cucCan + taeGut.calAnn + colLiv.mesUni + taeGut.colLiv + galGal.melGal + galGal.anaPla) - (2 * taeGut.galGal)) %>%
  mutate(neo_tip_loss = ctb(taeGut) + ctb(ficAlb) + ctb(pseHum) + ctb(corBra) + ctb(melUnd) + ctb(falPer) + ctb(picPub) + ctb(lepDis) + ctb(halLeu) + ctb(fulGla) + ctb(nipNip) + ctb(balReg) + ctb(chaVoc) + ctb(calAnn) + ctb(chaPel) + ctb(cucCan) + ctb(colLiv) + ctb(mesUni) + ctb(galGal) + ctb(melGal) +ctb(anaPla) + ctb(aptFor) + ctb(pygAde)) %>%
  mutate(tin_tip_loss = ctb(cryCin) + ctb(eudEle) + ctb(notPer) + ctb(tinGut)) %>%
  mutate(peng_pp_loss = (aptFor -  aptFor.pygAde ) + (pygAde - aptFor.pygAde) + (aptFor.pygAde - aptFor.fulGla)) %>%
  #now select columns to keep
  select(cnee, dataset, version, ends_with("_loss")) %>%
  #now compute total ratite losses with various assumptions
  mutate(floss_sp_pp = cd_pp_loss + rh_pp_loss + os_pp_loss + ki_pp_loss + mo_pp_loss) %>%
  mutate(floss_cl_pp = tto(cd_pp_loss) + tto(rh_pp_loss) + tto(os_pp_loss) + tto(ki_pp_loss) + tto(mo_pp_loss)) %>%
  mutate(floss_cl_pp_dollo = tto(cd_pp_loss + rh_pp_loss + ki_pp_loss) + tto(os_pp_loss) + tto(mo_pp_loss)) %>%
  mutate(vloss_pp = ti_pp_loss + it_pp_loss + neo_pp_loss + peng_pp_loss)

cnee_red_losses <- cnee_reduced %>%
  #this is the posterior estimated number of losses on each clade
  mutate(cd_pp_loss = (casCas - casCas.droNov) + (droNov - casCas.droNov) + (casCas.droNov - aptHaa.casCas)) %>%
  mutate(rh_pp_loss = (rheAme - rheAme.rhePen) + (rhePen - rheAme.rhePen) + (rheAme.rhePen - aptHaa.rheAme)) %>%
  mutate(os_pp_loss = strCam - aptHaa.strCam) %>%
  mutate(ki_pp_loss = (aptHaa - aptHaa.aptRow) + (aptOwe - aptHaa.aptOwe) + (aptRow - aptHaa.aptRow) + (aptHaa.aptRow - aptHaa.aptOwe) + (aptHaa.aptOwe - aptHaa.casCas)) %>%
  mutate(ti_pp_loss = (tinGut - cryCin.tinGut) + (cryCin - cryCin.tinGut) + (eudEle - eudEle.notPer) + (notPer - eudEle.notPer) + (eudEle.notPer - cryCin.eudEle) + (cryCin.tinGut - cryCin.eudEle)) %>%
  mutate(it_pp_loss = (cryCin.eudEle - aptHaa.cryCin) + (aptHaa.rheAme - aptHaa.cryCin) + (aptHaa.casCas - aptHaa.rheAme) + (aptHaa.cryCin - aptHaa.strCam)) %>%
  mutate(neo_pp_loss = (taeGut + ficAlb + pseHum + corBra + melUnd + falPer + picPub + lepDis + halLeu + fulGla + nipNip + balReg + chaVoc + calAnn + chaPel + cucCan +colLiv + mesUni + galGal + melGal + anaPla) - (taeGut.ficAlb + taeGut.pseHum + taeGut.corBra + taeGut.melUnd + taeGut.falPer + picPub.lepDis + picPub.halLeu + taeGut.picPub  + aptFor.nipNip + taeGut.aptFor + balReg.chaVoc + taeGut.balReg + calAnn.chaPel + calAnn.cucCan + taeGut.calAnn + colLiv.mesUni + taeGut.colLiv + galGal.melGal + galGal.anaPla) - (2 * taeGut.galGal)) %>%
  mutate(neo_tip_loss = ctb(taeGut) + ctb(ficAlb) + ctb(pseHum) + ctb(corBra) + ctb(melUnd) + ctb(falPer) + ctb(picPub) + ctb(lepDis) + ctb(halLeu) + ctb(fulGla) + ctb(nipNip) + ctb(balReg) + ctb(chaVoc) + ctb(calAnn) + ctb(chaPel) + ctb(cucCan) + ctb(colLiv) + ctb(mesUni) + ctb(galGal) + ctb(melGal) +ctb(anaPla) + ctb(aptFor) + ctb(pygAde)) %>%
  mutate(tin_tip_loss = ctb(cryCin) + ctb(eudEle) + ctb(notPer) + ctb(tinGut)) %>%
  mutate(peng_pp_loss = (aptFor -  aptFor.pygAde ) + (pygAde - aptFor.pygAde) + (aptFor.pygAde - aptFor.fulGla)) %>%
  #now select columns to keep
  select(cnee, dataset, version, ends_with("_loss")) %>%
  #now compute total ratite losses with various assumptions
  mutate(floss_sp_pp = cd_pp_loss + rh_pp_loss + os_pp_loss + ki_pp_loss) %>%
  mutate(floss_cl_pp = tto(cd_pp_loss) + tto(rh_pp_loss) + tto(os_pp_loss) + tto(ki_pp_loss)) %>%
  mutate(floss_cl_pp_dollo = tto(cd_pp_loss + rh_pp_loss + ki_pp_loss) + tto(os_pp_loss)) %>%
  mutate(vloss_pp = ti_pp_loss + it_pp_loss + neo_pp_loss + peng_pp_loss)

cnee_ext_losses <- cnee_ext %>%
  #this is the posterior estimated number of losses on each clade
  mutate(cd_pp_loss = (casCas - casCas.droNov) + (droNov - casCas.droNov) + (casCas.droNov - aptHaa.casCas)) %>%
  mutate(rh_pp_loss = (rheAme - rheAme.rhePen) + (rhePen - rheAme.rhePen) + (rheAme.rhePen - aptHaa.rheAme)) %>%
  mutate(os_pp_loss = strCam - aptHaa.strCam) %>%
  mutate(ki_pp_loss = (aptHaa - aptHaa.aptRow) + (aptOwe - aptHaa.aptOwe) + (aptRow - aptHaa.aptRow) + (aptHaa.aptRow - aptHaa.aptOwe) + (aptHaa.aptOwe - aptHaa.casCas)) %>%
  mutate(mo_pp_loss = anoDid - cryCin.anoDid) %>%
  mutate(ti_pp_loss = (tinGut - cryCin.tinGut) + (cryCin - cryCin.tinGut) + (eudEle - eudEle.notPer) + (notPer - eudEle.notPer) + (eudEle.notPer - cryCin.eudEle) + (cryCin.tinGut - cryCin.eudEle)) %>%
  mutate(it_pp_loss = (cryCin.eudEle - cryCin.anoDid) + (cryCin.anoDid - aptHaa.cryCin) + (aptHaa.rheAme - aptHaa.cryCin) + (aptHaa.casCas - aptHaa.rheAme) + (aptHaa.cryCin - aptHaa.strCam)) %>%
  mutate(neo_pp_loss = (taeGut + ficAlb + pseHum + corBra + melUnd + falPer + picPub + lepDis + halLeu + fulGla + nipNip + balReg + chaVoc + calAnn + chaPel + cucCan +colLiv + mesUni + galGal + melGal + anaPla + nanAur + nanBra + uriPel) - (taeGut.ficAlb + taeGut.pseHum + taeGut.corBra + taeGut.melUnd + taeGut.falPer + picPub.lepDis + picPub.halLeu + taeGut.picPub  + taeGut.aptFor + balReg.chaVoc + taeGut.balReg + calAnn.chaPel + calAnn.cucCan + taeGut.calAnn + colLiv.mesUni + taeGut.colLiv + galGal.melGal + galGal.anaPla + nanAur.nanBra + nanAur.uriPel + nanAur.nipNip + aptFor.nanAur) - (2 * taeGut.galGal)) %>%
  mutate(neo_tip_loss = ctb(taeGut) + ctb(ficAlb) + ctb(pseHum) + ctb(corBra) + ctb(melUnd) + ctb(falPer) + ctb(picPub) + ctb(lepDis) + ctb(halLeu) + ctb(fulGla) + ctb(nipNip) + ctb(balReg) + ctb(chaVoc) + ctb(calAnn) + ctb(chaPel) + ctb(cucCan) + ctb(colLiv) + ctb(mesUni) + ctb(galGal) + ctb(melGal) +ctb(anaPla) + ctb(aptFor) + ctb(pygAde) + ctb(nanAur) + ctb(nanBra) + ctb(uriPel)) %>%
  mutate(tin_tip_loss = ctb(cryCin) + ctb(eudEle) + ctb(notPer) + ctb(tinGut)) %>%
  mutate(peng_pp_loss = (aptFor -  aptFor.pygAde ) + (pygAde - aptFor.pygAde) + (aptFor.pygAde -  aptFor.fulGla)) %>%
  mutate(gc_pp_loss = nanHar - nanAur.nanHar) %>%
  #now select columns to keep
  select(cnee, dataset, version, ends_with("_loss")) %>%
  #now compute total ratite losses with various assumptions
  mutate(floss_sp_pp = cd_pp_loss + rh_pp_loss + os_pp_loss + ki_pp_loss + mo_pp_loss + gc_pp_loss) %>%
  mutate(floss_cl_pp = tto(cd_pp_loss) + tto(rh_pp_loss) + tto(os_pp_loss) + tto(ki_pp_loss) + tto(mo_pp_loss) + tto(gc_pp_loss)) %>%
  mutate(floss_cl_pp_dollo = tto(cd_pp_loss + rh_pp_loss + ki_pp_loss) + tto(os_pp_loss) + tto(mo_pp_loss) + tto(gc_pp_loss)) %>%
  mutate(vloss_pp = ti_pp_loss + it_pp_loss + neo_pp_loss + peng_pp_loss)

#write out final datasets with everything
cnee_orig_final <- inner_join(cnee_orig, cnee_orig_losses, by=c("cnee" = "cnee", "version" = "version", "dataset" = "dataset"))
cnee_red_final <- inner_join(cnee_reduced, cnee_red_losses, by=c("cnee" = "cnee", "version" = "version", "dataset" = "dataset"))
cnee_ext_final <- inner_join(cnee_ext, cnee_ext_losses, by=c("cnee" = "cnee", "version" = "version", "dataset" = "dataset"))

#write_tsv(cnee_orig_final, path="/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/cnees_original.tsv")
#write_tsv(cnee_ext_final, path="/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/cnees_extended.tsv")
#write_tsv(cnee_red_final, path="/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/cnees_reduced.tsv")

## BASIC QC ##

setwd("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/04_qc")

#read in original runs
final_v1 <- read_tsv("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/cnees.tsv.gz")

#bf2 vs nonratite loss
cnee_ext_final %>% filter(logBF1 >= 10) %>% 
  mutate(tip_loss = ifelse(neo_tip_loss+tin_tip_loss > 10, 10, neo_tip_loss+tin_tip_loss)) %>%
  ggplot(aes(as.factor(tip_loss), logBF2, fill=version)) + 
  geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim=c(-60,20)) 
ggsave("bf2_vs_tiploss_extended.pdf")

cnee_red_final %>% filter(logBF1 >= 10) %>% 
  mutate(tip_loss = ifelse(neo_tip_loss+tin_tip_loss > 10, 10, neo_tip_loss+tin_tip_loss)) %>%
  ggplot(aes(as.factor(tip_loss), logBF2, fill=version)) + 
  geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim=c(-60,20)) 
ggsave("bf2_vs_tiploss_reduced.pdf")

cnee_orig_final %>% filter(logBF1 >= 10) %>% 
  mutate(tip_loss = ifelse(neo_tip_loss+tin_tip_loss > 10, 10, neo_tip_loss+tin_tip_loss)) %>%
  ggplot(aes(as.factor(tip_loss), logBF2, fill=version)) + 
  geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim=c(-60,20)) 
ggsave("bf2_vs_tiploss_orig.pdf")

#prob loss in ostrich between original runs and new runs

cnee_red_final %>% select(cnee, version, strCam) %>% full_join(final_v1, by=c("cnee" = "cnee")) %>% 
  filter(version == "orig") %>%
  ggplot(aes(strCam, strCam.prob)) + geom_point(alpha=0.1) + geom_abline(col="red", size=1, linetype="dashed")

cnee_red_final %>% select(cnee, version, strCam) %>% full_join(final_v1, by=c("cnee" = "cnee")) %>% 
  filter(version == "gain") %>%
  ggplot(aes(strCam, strCam.prob)) + geom_point(alpha=0.1) + geom_abline(col="red", size=1, linetype="dashed")


#distribution of # of losses for different datasets

cnee_orig_final %>% filter(logBF1 >= 10, logBF2 >= 1) %>% 
  ggplot(aes(floss_cl_pp, color=version)) +
  geom_density()
  
cnee_orig_final %>% filter(logBF1 >= 10, logBF2 >= 1, tin_tip_loss == 0) %>% 
  ggplot(aes(floss_cl_pp, color=version)) +
  geom_density()

cnee_orig_final %>% filter(logBF1 >= 10, logBF2 >= 1, tin_tip_loss == 0, neo_tip_loss == 0) %>% 
  ggplot(aes(floss_cl_pp, color=version)) +
  geom_density()

final_v1 %>% filter(bf1 >= 10, bf2 >= 1) %>% 
  ggplot(aes(ratite_loss_cons.prob)) +
  geom_density()

cnee_orig_final %>% filter(logBF1 >= 10, logBF2 >= 1) %>% group_by(version) %>%
  mutate(conv = ifelse(floss_cl_pp >= 2, "conv", "sing")) %>% count(conv) %>% spread(conv, n) %>%
  mutate(frac = conv / (conv+sing))

cnee_orig_final %>% filter(logBF1 >= 10, logBF2 >= 1, tin_tip_loss == 0) %>% group_by(version) %>%
  mutate(conv = ifelse(floss_cl_pp >= 2, "conv", "sing")) %>% count(conv) %>% spread(conv, n) %>%
  mutate(frac = conv / (conv+sing))

cnee_orig_final %>% filter(logBF1 >= 10, logBF2 >= 1, tin_tip_loss == 0, neo_tip_loss == 0) %>% 
  group_by(version) %>%
  mutate(conv = ifelse(floss_cl_pp >= 2, "conv", "sing")) %>% count(conv) %>% spread(conv, n) %>%
  mutate(frac = conv / (conv+sing))

final_v1 %>% filter(bf1 >= 10, bf2 >= 1) %>%
  mutate(conv = ifelse(ratite_loss_cons.prob >= 2, "conv", "sing")) %>% count(conv) %>% spread(conv, n) %>%
  mutate(frac = conv / (conv+sing))


#let's just compare final_v1 loss prob to floss_cl_pp for ratite specific

cnee_orig_final %>% filter(version=="orig") %>%
  select(cnee, logBF1, logBF2, floss_cl_pp) %>% inner_join(final_v1) %>%
  filter(logBF1 >= 10 & logBF2 >= 1 | bf1 >= 10 & bf2 >= 1) %>%
  ggplot(aes(floss_cl_pp, ratite_loss_cons.prob)) + geom_point(alpha=0.1) +
  geom_abline(col="red", size=1, linetype="dashed")
  

cnee_orig_final %>% filter(version=="gap") %>%
  select(cnee, logBF1, logBF2, floss_cl_pp) %>% inner_join(final_v1) %>%
  filter(logBF1 >= 10 & logBF2 >= 1 | bf1 >= 10 & bf2 >= 1) %>%
  ggplot(aes(floss_cl_pp, ratite_loss_cons.prob)) + geom_point(alpha=0.1) +
  xlab("Ratite Losses, New Run") + 
  ylab("Ratite Losses, Old Run") +
  geom_abline(col="red", size=1, linetype="dashed")


cnee_orig_final %>% filter(version=="gap") %>%
  select(cnee, logBF1, logBF2, floss_cl_pp, strCam, anoDid, rheAme, rhePen, droNov, casCas) %>% inner_join(final_v1) %>%
  filter(logBF1 >= 10 & logBF2 >= 1 | bf1 >= 10 & bf2 >= 1) %>%
  ggplot(aes(droNov, droNov.prob)) + geom_point(alpha=0.1) +
  geom_abline(col="red", size=1, linetype="dashed")


#convergence by tip losses only, comparing ratites and random

cnee_orig_final %>% filter(version == "orig") %>% 
  filter(tin_tip_loss == 0, neo_tip_loss == 0) %>%
  select(cnee, strCam, anoDid, casCas) %>%
  mutate(conv_ct = ctb(strCam) + ctb(anoDid) + ctb(casCas)) %>%
  count(conv_ct)
  
cnee_orig_final %>% filter(version == "orig") %>% 
  select(cnee, taeGut:strCam, tin_tip_loss, neo_tip_loss) %>%
  mutate(total_tip_loss = ctb(strCam) + ctb(anoDid) + ctb(casCas) + ctb(droNov) + ctb(aptHaa) + ctb(aptOwe) + + ctb(aptRow) + ctb(rheAme) + ctb(rhePen) + tin_tip_loss + neo_tip_loss) %>%
  mutate(conv_ct = ctb(balReg) + ctb(calAnn) + ctb(pygAde)) %>%
  filter(conv_ct == total_tip_loss) %>%
  count(conv_ct)


cnee_orig_final %>% filter(version=="gap") %>%
  select(cnee, logBF1, logBF2, floss_cl_pp, strCam, anoDid, rheAme, rhePen, droNov, casCas) %>% inner_join(final_v1) %>%
  filter(logBF1 >= 10 & logBF2 >= 1 | bf1 >= 10 & bf2 >= 1) %>%
  filter(droNov < 0.25, droNov.prob > 0.99) %>% select(cnee, logBF1, logBF2, bf1, bf2)


cnee_orig_final %>% filter(version=="gap" | version=="orig", logBF1 >= 10, logBF2 >= 1)%>%
  select(cnee, version, floss_cl_pp) %>%
  spread(version, floss_cl_pp) %>%
  ggplot(aes(gap, orig)) + geom_point(alpha=0.1) +
  geom_abline(col="red", size=1, linetype="dashed")


#basic QC looks okay

