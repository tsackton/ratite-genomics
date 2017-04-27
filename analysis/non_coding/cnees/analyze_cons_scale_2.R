#confirm Z matrix results
setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/")
library(dplyr)
library(tidyr)
library(data.table)

#read in original analysis
bayes_res <- read.table("all_bayes_results.tsv", stringsAsFactors = F) %>% tibble::rownames_to_column(var="cnee") %>% tbl_df

bayes_lik <- read.table("/Volumes/LaCie/Projects/Current/ratites/final/cons_scale/version2/subseq_cons_ALL_elem_lik.txt", stringsAsFactors = F, header=F) %>% tbl_df %>% rename(num = V1,	cnee = V2,	ll_null = V3,	ll_resZ = V4,	ll_full = V5, bf = V6, ll_max1 = V7, ll_max2 = V8) %>% select(-num) %>% inner_join(., bayes_res)

bayes_long <- bayes_lik %>% gather(species, accel, taeGut:strCam)

get_clade <- function(x, type) {
  if (x %in% c("rheAme", "rhePen")) {
    if (type == "bg") { return(c("rhea")) }
    return(c("ekcr"))
  }
  if (x %in% c("strCam")) {
    return(c("ostrich"))
  }
  if (x %in% c("anoDid")) {
    return(c("moa"))
  }
  if (x %in% c("aptRow", "aptHaa", "aptOwe")) {
    if (type == "bg") { return(c("kiwi")) }
    return(c("ekcr"))
  }
  if (x %in% c("droNov", "casCas")) {
    if (type == "bg") { return(c("emucas"))}
    return(c("ekcr"))
  }
  return(c("non_ratite"))
}

tipnames <- distinct(bayes_long, species) %>% group_by(species) %>% mutate(clade_key_bg = get_clade(species, "bg"), clade_key_par = get_clade(species, "par"))

bayes_bg_clade <- bayes_long %>% inner_join(., tipnames) %>% group_by(cnee, clade_key_bg) %>% mutate(clade_accel = max(accel)) %>% ungroup %>% select(cnee, bf, clade_key_bg, clade_accel) %>% distinct %>% spread(clade_key_bg, clade_accel) %>% mutate(ratite_ct = emucas + kiwi + moa + ostrich + rhea)
bayes_bg_clade %>% filter(non_ratite == 0) %>% with(., table(bf > 1, ratite_ct))

bayes_bg_clade %>% filter(non_ratite == 0) %>% group_by(ratite_ct) %>% summarize(median=median(bf), total=n(), count=sum(bf>0), count_strong=sum(bf>3)) %>% mutate(frac_pos = count/total, frac_strong = count_strong/total)

bayes_bg_clade %>% arrange(desc(bf)) %>% head(n=500) %>% with(., table(ratite_ct, non_ratite))

bayes_bg_clade %>% filter(non_ratite == 0, ratite_ct >= 4) %>% summarize(median=median(bf))
bayes_bg_clade %>% filter(non_ratite == 0, ratite_ct >= 3) %>% summarize(median=median(bf))


#parsimony
bayes_par_clade <- bayes_long %>% inner_join(., tipnames) %>% group_by(cnee, clade_key_par) %>% mutate(clade_accel = max(accel)) %>% ungroup %>% select(cnee, bf, clade_key_par, clade_accel) %>% distinct %>% spread(clade_key_par, clade_accel) %>% mutate(ratite_ct = ekcr + moa + ostrich)
bayes_par_clade %>% filter(non_ratite == 0) %>% with(., table(bf > 1, ratite_ct))

bayes_par_clade %>% filter(non_ratite == 0) %>% group_by(ratite_ct) %>% summarize(median=median(bf), total=n(), count=sum(bf>0), count_strong=sum(bf>3)) %>% mutate(frac_pos = count/total, frac_strong = count_strong/total)

bayes_par_clade %>% arrange(desc(bf)) %>% head(n=500) %>% with(., table(ratite_ct, non_ratite))

bayes_par_clade %>% filter(non_ratite == 0, ratite_ct >= 3) %>% summarize(median=median(bf))


#COMPARE TO PHAST

#read in phast
load("/Volumes/LaCie/Projects/Current/ratites/final/accelTests/all.accel.processed.Rdata")

#get just the "withMoa" runtype and neut_ver3
phast<-all.accel %>% filter(runtype=="withMoa", model=="neut_ver3")
phast<-phast %>% group_by(group) %>% mutate(qval=p.adjust(pval, method="fdr")) %>% ungroup
phast_sum<-phast %>% mutate(sig = as.numeric(qval < 0.05)) %>% select(name, group, sig) %>% distinct %>%
  spread(group, sig)

phast_sum <- phast_sum %>% mutate(ratite_clade = anoDid + strCam + Kiwi + Casuar + Rhea) %>% select(cnee = name, phast_ct = ratite_clade)

comp<-inner_join(bayes_bg_clade, phast_sum)

comp %>% filter(non_ratite == 0) %>% with(., table(ratite_ct, phast_ct >= 1))

comp %>% filter(bf > 1, non_ratite == 0, phast_ct >=3, ratite_ct >= 3, moa == 1 | ostrich == 1) %>% inner_join(.,phast, by=c("cnee" = "name")) %>% filter(group == "Tinamou" | group == "BasalPaleo", qval > 0.25) %>% distinct(cnee)
