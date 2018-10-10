## CODE TO PARSE PHYLOACC OUTPUT ##
## UPDATED OCT 2018 FOR MANUSCRIPT REVISIONS ##

library(tidyverse)
setwd("~/Projects/birds/ratite_compgen/data/2018-09-FINAL/CNEE/cnee_results/")

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

input_lik <- function(file, dataset) {
  read_delim(file=file, delim="\t") %>% 
    mutate(bf1 = loglik_RES - loglik_NUll, bf2 = loglik_RES - loglik_all) %>%
    mutate(dataset = dataset) %>%
    select(dataset, cnee = ID, bf1, bf2, ll_null = loglik_NUll, ll_target = loglik_RES, ll_full = loglik_all)
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
  posteriors[[dataset]] %>% 
  select(cnee, n_rate, c_rate, g_rate, l_rate, l2_rate, ends_with("_3")) %>%
  rename_at(vars(ends_with("_3")), gsub, pattern="_3", replacement="", fixed=TRUE) %>% 
  inner_join(likelihoods[[dataset]], ., by=c("cnee" = "cnee")) %>% 
  select(-dataset) %>%
  inner_join(postmat[[dataset]], by=c("cnee" = "cnee"), suffix = c("", ".max")) %>%
  select(dataset, cnee, everything())
}

## END FUNCTIONS ##

likelihoods <- list()
posteriors <- list()
postmat <- list()

#read in all files; not automated due to variation in file names, etc
likelihoods$original <- input_lik(file="/Users/tim/Projects/birds/ratite_compgen/data/2018-09-FINAL/CNEE/cnee_results/phyloAcc_output/original_set/Combined_elem_lik2-1.txt", "original")
likelihoods$version1 <- input_lik(file="/Users/tim/Projects/birds/ratite_compgen/data/2018-09-FINAL/CNEE/cnee_results/phyloAcc_output/original_set/phyloAcc_v1/element_likelihood.txt", "version1")
likelihoods$extended <- input_lik(file="/Users/tim/Projects/birds/ratite_compgen/data/2018-09-FINAL/CNEE/cnee_results/phyloAcc_output/extended_set/Combined_elem_lik-E.txt", "extended")
likelihoods$reduced <- input_lik(file="/Users/tim/Projects/birds/ratite_compgen/data/2018-09-FINAL/CNEE/cnee_results/phyloAcc_output/reduced_set/Combined_elem_lik-R.txt", "reduced")


postmat$original <- input_states(file="/Users/tim/Projects/birds/ratite_compgen/data/2018-09-FINAL/CNEE/cnee_results/phyloAcc_output/original_set/Combined_max_Z-1.txt", "original", "maxpost")
postmat$version1 <- read_tsv(file="/Users/tim/Projects/birds/ratite_compgen/data/2018-09-FINAL/CNEE/cnee_results/phyloAcc_output/original_set/phyloAcc_v1/max_Z.txt") %>% mutate(cnee = likelihoods$version1$cnee)
postmat$extended <- input_states(file="/Users/tim/Projects/birds/ratite_compgen/data/2018-09-FINAL/CNEE/cnee_results/phyloAcc_output/extended_set/Combined_max_Z-E.txt", "extended", "maxpost")
postmat$reduced <- input_states(file="/Users/tim/Projects/birds/ratite_compgen/data/2018-09-FINAL/CNEE/cnee_results/phyloAcc_output/reduced_set/Combined_max_Z-R.txt", "reduced", "maxpost")

posteriors$original <- input_states(file="/Users/tim/Projects/birds/ratite_compgen/data/2018-09-FINAL/CNEE/cnee_results/phyloAcc_output/original_set/Combined_post_Z_1-1.txt", "original", "postprob")
posteriors$version1 <- read_tsv(file="/Users/tim/Projects/birds/ratite_compgen/data/2018-09-FINAL/CNEE/cnee_results/phyloAcc_output/original_set/phyloAcc_v1/post_Z.txt") %>% mutate(cnee = likelihoods$version1$cnee)
posteriors$extended <- input_states(file="/Users/tim/Projects/birds/ratite_compgen/data/2018-09-FINAL/CNEE/cnee_results/phyloAcc_output/extended_set/Combined_post_Z_1-E.txt", "extended", "postprob")
posteriors$reduced <- input_states(file="/Users/tim/Projects/birds/ratite_compgen/data/2018-09-FINAL/CNEE/cnee_results/phyloAcc_output/reduced_set/Combined_post_Z_1-R.txt", "reduced", "postprob")

#now merge each dataset to make one cnee set for each phyloAcc run, containing the CNEE id, basic parameters, post prob acceleration, max a posteriori estimated state

cnee_orig <- merge_phyloAcc("original")
cnee_reduced <- merge_phyloAcc("reduced")
cnee_ext <- merge_phyloAcc("extended")
cnee_v1 <- posteriors$version1 %>% 
    select(cnee, n_rate, c_rate, ends_with("_3")) %>%
    setNames(c("cnee", "n_rate", "c_rate", colnames(postmat$version1)[1:length(colnames(postmat$version1))-1])) %>% 
    inner_join(likelihoods$version1, ., by=c("cnee" = "cnee")) %>% 
    inner_join(postmat$version1, by=c("cnee" = "cnee"), suffix = c("", ".max")) %>%
    select(dataset, cnee, everything())


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
  mutate(nv_pp_loss = (taeGut + ficAlb + pseHum + corBra + melUnd + falPer + picPub + lepDis + halLeu + fulGla + nipNip + balReg + chaVoc + calAnn + chaPel + cucCan +colLiv + mesUni + galGal + melGal + anaPla) - (taeGut.ficAlb + taeGut.pseHum + taeGut.corBra + taeGut.melUnd + taeGut.falPer + picPub.lepDis + picPub.halLeu + taeGut.picPub  + aptFor.nipNip + taeGut.aptFor + balReg.chaVoc + taeGut.balReg + calAnn.chaPel + calAnn.cucCan + taeGut.calAnn + colLiv.mesUni + taeGut.colLiv + galGal.melGal + galGal.anaPla) - (2 * taeGut.galGal)) %>%
  mutate(pg_pp_loss = (aptFor -  aptFor.pygAde ) + (pygAde - aptFor.pygAde) + (aptFor.pygAde - aptFor.fulGla)) %>%
  #now select columns to keep
  select(cnee, ends_with("_loss")) %>%
  #now compute total ratite losses with various assumptions
  mutate(floss_sp_pp = cd_pp_loss + rh_pp_loss + os_pp_loss + ki_pp_loss + mo_pp_loss) %>%
  mutate(floss_cl_pp = tto(cd_pp_loss) + tto(rh_pp_loss) + tto(os_pp_loss) + tto(ki_pp_loss) + tto(mo_pp_loss)) %>%
  mutate(floss_cl_pp_dollo = tto(cd_pp_loss + rh_pp_loss + ki_pp_loss) + tto(os_pp_loss) + tto(mo_pp_loss)) %>%
  mutate(vloss_pp = ti_pp_loss + it_pp_loss + nv_pp_loss + pg_pp_loss)

cnee_v1_losses <- cnee_v1 %>%
  #this is the posterior estimated number of losses on each clade
  mutate(cd_pp_loss = (casCas - casCas.droNov) + (droNov - casCas.droNov) + (casCas.droNov - aptHaa.casCas)) %>%
  mutate(rh_pp_loss = (rheAme - rheAme.rhePen) + (rhePen - rheAme.rhePen) + (rheAme.rhePen - aptHaa.rheAme)) %>%
  mutate(os_pp_loss = strCam - aptHaa.strCam) %>%
  mutate(ki_pp_loss = (aptHaa - aptHaa.aptRow) + (aptOwe - aptHaa.aptOwe) + (aptRow - aptHaa.aptRow) + (aptHaa.aptRow - aptHaa.aptOwe) + (aptHaa.aptOwe - aptHaa.casCas)) %>%
  mutate(mo_pp_loss = anoDid - cryCin.anoDid) %>%
  mutate(ti_pp_loss = (tinGut - cryCin.tinGut) + (cryCin - cryCin.tinGut) + (eudEle - eudEle.notPer) + (notPer - eudEle.notPer) + (eudEle.notPer - cryCin.eudEle) + (cryCin.tinGut - cryCin.eudEle)) %>%
  mutate(it_pp_loss = (cryCin.eudEle - cryCin.anoDid) + (cryCin.anoDid - aptHaa.cryCin) + (aptHaa.rheAme - aptHaa.cryCin) + (aptHaa.casCas - aptHaa.rheAme) + (aptHaa.cryCin - aptHaa.strCam)) %>%
  mutate(nv_pp_loss = (taeGut + ficAlb + pseHum + corBra + melUnd + falPer + picPub + lepDis + halLeu + fulGla + nipNip + balReg + chaVoc + calAnn + chaPel + cucCan +colLiv + mesUni + galGal + melGal + anaPla) - (taeGut.ficAlb + taeGut.pseHum + taeGut.corBra + taeGut.melUnd + taeGut.falPer + picPub.lepDis + picPub.halLeu + taeGut.picPub  + aptFor.nipNip + taeGut.aptFor + balReg.chaVoc + taeGut.balReg + calAnn.chaPel + calAnn.cucCan + taeGut.calAnn + colLiv.mesUni + taeGut.colLiv + galGal.melGal + galGal.anaPla) - (2 * taeGut.galGal)) %>%
  mutate(pg_pp_loss = (aptFor - aptFor.pygAde ) + (pygAde - aptFor.pygAde) + (aptFor.pygAde - aptFor.fulGla)) %>%
  #now select columns to keep
  select(cnee, ends_with("_loss")) %>%
  #now compute total ratite losses with various assumptions
  mutate(floss_sp_pp = cd_pp_loss + rh_pp_loss + os_pp_loss + ki_pp_loss + mo_pp_loss) %>%
  mutate(floss_cl_pp = tto(cd_pp_loss) + tto(rh_pp_loss) + tto(os_pp_loss) + tto(ki_pp_loss) + tto(mo_pp_loss)) %>%
  mutate(floss_cl_pp_dollo = tto(cd_pp_loss + rh_pp_loss + ki_pp_loss) + tto(os_pp_loss) + tto(mo_pp_loss)) %>%
  mutate(vloss_pp = ti_pp_loss + it_pp_loss + nv_pp_loss + pg_pp_loss)

cnee_red_losses <- cnee_reduced %>%
  #this is the posterior estimated number of losses on each clade
  mutate(cd_pp_loss = (casCas - casCas.droNov) + (droNov - casCas.droNov) + (casCas.droNov - aptHaa.casCas)) %>%
  mutate(rh_pp_loss = (rheAme - rheAme.rhePen) + (rhePen - rheAme.rhePen) + (rheAme.rhePen - aptHaa.rheAme)) %>%
  mutate(os_pp_loss = strCam - aptHaa.strCam) %>%
  mutate(ki_pp_loss = (aptHaa - aptHaa.aptRow) + (aptOwe - aptHaa.aptOwe) + (aptRow - aptHaa.aptRow) + (aptHaa.aptRow - aptHaa.aptOwe) + (aptHaa.aptOwe - aptHaa.casCas)) %>%
  mutate(ti_pp_loss = (tinGut - cryCin.tinGut) + (cryCin - cryCin.tinGut) + (eudEle - eudEle.notPer) + (notPer - eudEle.notPer) + (eudEle.notPer - cryCin.eudEle) + (cryCin.tinGut - cryCin.eudEle)) %>%
  mutate(it_pp_loss = (cryCin.eudEle - aptHaa.cryCin) + (aptHaa.rheAme - aptHaa.cryCin) + (aptHaa.casCas - aptHaa.rheAme) + (aptHaa.cryCin - aptHaa.strCam)) %>%
  mutate(nv_pp_loss = (taeGut + ficAlb + pseHum + corBra + melUnd + falPer + picPub + lepDis + halLeu + fulGla + nipNip + balReg + chaVoc + calAnn + chaPel + cucCan +colLiv + mesUni + galGal + melGal + anaPla) - (taeGut.ficAlb + taeGut.pseHum + taeGut.corBra + taeGut.melUnd + taeGut.falPer + picPub.lepDis + picPub.halLeu + taeGut.picPub  + aptFor.nipNip + taeGut.aptFor + balReg.chaVoc + taeGut.balReg + calAnn.chaPel + calAnn.cucCan + taeGut.calAnn + colLiv.mesUni + taeGut.colLiv + galGal.melGal + galGal.anaPla) - (2 * taeGut.galGal)) %>%
  mutate(pg_pp_loss = (aptFor -  aptFor.pygAde ) + (pygAde - aptFor.pygAde) + (aptFor.pygAde - aptFor.fulGla)) %>%
  #now select columns to keep
  select(cnee, ends_with("_loss")) %>%
  #now compute total ratite losses with various assumptions
  mutate(floss_sp_pp = cd_pp_loss + rh_pp_loss + os_pp_loss + ki_pp_loss) %>%
  mutate(floss_cl_pp = tto(cd_pp_loss) + tto(rh_pp_loss) + tto(os_pp_loss) + tto(ki_pp_loss)) %>%
  mutate(floss_cl_pp_dollo = tto(cd_pp_loss + rh_pp_loss + ki_pp_loss) + tto(os_pp_loss)) %>%
  mutate(vloss_pp = ti_pp_loss + it_pp_loss + nv_pp_loss + pg_pp_loss)

cnee_ext_losses <- cnee_ext %>%
  #this is the posterior estimated number of losses on each clade
  mutate(cd_pp_loss = (casCas - casCas.droNov) + (droNov - casCas.droNov) + (casCas.droNov - aptHaa.casCas)) %>%
  mutate(rh_pp_loss = (rheAme - rheAme.rhePen) + (rhePen - rheAme.rhePen) + (rheAme.rhePen - aptHaa.rheAme)) %>%
  mutate(os_pp_loss = strCam - aptHaa.strCam) %>%
  mutate(ki_pp_loss = (aptHaa - aptHaa.aptRow) + (aptOwe - aptHaa.aptOwe) + (aptRow - aptHaa.aptRow) + (aptHaa.aptRow - aptHaa.aptOwe) + (aptHaa.aptOwe - aptHaa.casCas)) %>%
  mutate(mo_pp_loss = anoDid - cryCin.anoDid) %>%
  mutate(ti_pp_loss = (tinGut - cryCin.tinGut) + (cryCin - cryCin.tinGut) + (eudEle - eudEle.notPer) + (notPer - eudEle.notPer) + (eudEle.notPer - cryCin.eudEle) + (cryCin.tinGut - cryCin.eudEle)) %>%
  mutate(it_pp_loss = (cryCin.eudEle - cryCin.anoDid) + (cryCin.anoDid - aptHaa.cryCin) + (aptHaa.rheAme - aptHaa.cryCin) + (aptHaa.casCas - aptHaa.rheAme) + (aptHaa.cryCin - aptHaa.strCam)) %>%
  mutate(nv_pp_loss = (taeGut + ficAlb + pseHum + corBra + melUnd + falPer + picPub + lepDis + halLeu + fulGla + nipNip + balReg + chaVoc + calAnn + chaPel + cucCan +colLiv + mesUni + galGal + melGal + anaPla + nanAur + nanBra + uriPel) - (taeGut.ficAlb + taeGut.pseHum + taeGut.corBra + taeGut.melUnd + taeGut.falPer + picPub.lepDis + picPub.halLeu + taeGut.picPub  + taeGut.aptFor + balReg.chaVoc + taeGut.balReg + calAnn.chaPel + calAnn.cucCan + taeGut.calAnn + colLiv.mesUni + taeGut.colLiv + galGal.melGal + galGal.anaPla + nanAur.nanBra + nanAur.uriPel + nanAur.nipNip + aptFor.nanAur) - (2 * taeGut.galGal)) %>%
  mutate(pg_pp_loss = (aptFor -  aptFor.pygAde ) + (pygAde - aptFor.pygAde) + (aptFor.pygAde -  aptFor.fulGla)) %>%
  mutate(gc_pp_loss = nanHar - nanAur.nanHar) %>%
  #now select columns to keep
  select(cnee, ends_with("_loss")) %>%
  #now compute total ratite losses with various assumptions
  mutate(floss_sp_pp = cd_pp_loss + rh_pp_loss + os_pp_loss + ki_pp_loss + mo_pp_loss + gc_pp_loss) %>%
  mutate(floss_cl_pp = tto(cd_pp_loss) + tto(rh_pp_loss) + tto(os_pp_loss) + tto(ki_pp_loss) + tto(mo_pp_loss) + tto(gc_pp_loss)) %>%
  mutate(floss_cl_pp_dollo = tto(cd_pp_loss + rh_pp_loss + ki_pp_loss) + tto(os_pp_loss) + tto(mo_pp_loss) + tto(gc_pp_loss)) %>%
  mutate(vloss_pp = ti_pp_loss + it_pp_loss + nv_pp_loss + pg_pp_loss)

#write out final datasets with everything
cnee_orig_final <- inner_join(cnee_orig, cnee_orig_losses, by=c("cnee" = "cnee"))
cnee_red_final <- inner_join(cnee_reduced, cnee_red_losses, by=c("cnee" = "cnee"))
cnee_ext_final <- inner_join(cnee_ext, cnee_ext_losses, by=c("cnee" = "cnee"))
cnee_v1_final <- inner_join(cnee_v1, cnee_v1_losses, by=c("cnee" = "cnee"))

write_tsv(cnee_orig_final, path="cnees_original.tsv")
write_tsv(cnee_ext_final, path="cnees_extended.tsv")
write_tsv(cnee_red_final, path="cnees_reduced.tsv")
write_tsv(cnee_v1_final, path="cnees_version1.tsv")
