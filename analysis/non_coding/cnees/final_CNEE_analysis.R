#new CNEE analysis - final version
#May 2017

library(tidyr)
library(dplyr)
library(ggplot2)

#read in data and clean up -- phyloAcc runs
setwd("/Volumes/LaCie/Projects/Current/ratites/final/phyloAcc_May2017/")
lik<-read.table("Combined_elem_lik_06.txt", header=T, stringsAsFactors = F) %>% tbl_df
zpost<-read.table("Combined_post_Z_06.txt", header=T, stringsAsFactors = F) %>% tbl_df
zmat<-read.table("Combined_max_Z_06.txt", header=T, stringsAsFactors = F) %>% tbl_df
zpost_clean <- zpost %>% 

#read in data and clean up -- phyloP runs
setwd("/Volumes/LaCie/Projects/Current/ratites/final/accelTests/withMoa/")
allRatite <- read.table("allRatite.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df
anoDid <- read.table("anoDid.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
casuar <- read.table("Casuar.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
kiwi <- read.table("Kiwi.out_neut_ver3.results", header =T , stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
rhea <- read.table("Rhea.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))
strCam <- read.table("strCam.out_neut_ver3.results", header = T, stringsAsFactors = F) %>% tbl_df %>% mutate(qval = p.adjust(pval, method="fdr"))

#merge and clean up
cnee.res <- lik %>% full_join(., allRatite, by = c("ID" = "name")) %>% select(ID, loglik_null = loglik_NUll, loglik_ratite = loglik_RES, loglik_full = loglik_all, phyloP_rate = null_scale, phyloP_altrate = alt_scale, phyloP_accrate = alt_subscale, phyloP_lr = lnlratio, phyloP_pval = pval) %>% mutate(phyloP_qval = p.adjust(phyloP_pval, method="fdr"), bf_null = loglik_ratite - loglik_null, bf_full = loglik_ratite - loglik_full)
cnee.res <- zpost %>% select(phyloAcc_accrate = n_rate, phyloAcc_rate = c_rate) %>% mutate(ID = lik$ID) %>% full_join(cnee.res, ., by=c("ID" = "ID"))

#support is highly correlated
with(cnee.res, boxplot(phyloP_lr ~ cut(bf_null, breaks=c(-1000,-10,0,1,5,10,20,1000)), outline=F))
cnee.res %>% ggplot(aes(x=phyloP_lr, y=bf_null)) + stat_binhex(bins=100)
fisher.test(table(cnee.res$phyloP_qval < 0.01, cnee.res$bf_null > 0))
fisher.test(table(cnee.res$phyloP_qval < 0.01, cnee.res$bf_null > 5))
fisher.test(table(cnee.res$phyloP_qval < 0.01, cnee.res$bf_null > 15))
fisher.test(table(cnee.res$phyloP_qval < 0.05, cnee.res$bf_null > 0))
fisher.test(table(cnee.res$phyloP_qval < 0.05, cnee.res$bf_null > 5))
fisher.test(table(cnee.res$phyloP_qval < 0.05, cnee.res$bf_null > 15))
fisher.test(table(cnee.res$phyloP_qval < 0.1, cnee.res$bf_null > 0))
fisher.test(table(cnee.res$phyloP_qval < 0.1, cnee.res$bf_null > 5))
fisher.test(table(cnee.res$phyloP_qval < 0.1, cnee.res$bf_null > 15))

#rates are correlated
cnee.res %>% ggplot(aes(x=phyloP_altrate, y=phyloAcc_rate)) + stat_binhex(bins=100) + scale_x_log10() + scale_y_log10()

cnee.res %>% ggplot(aes(x=phyloP_accrate, y=phyloAcc_accrate)) + stat_binhex(bins=100) + scale_x_log10() + scale_y_log10()

#classify as ratite-accelerated
cnee.res <- cnee.res %>% mutate(phyloAcc_ra = ifelse(bf_null > 5, TRUE, FALSE), phyloP_ra = ifelse(phyloP_qval < 0.01, TRUE, FALSE))
boxplot(cnee.res$phyloP_lr ~ cnee.res$phyloAcc_ra, outline=F)
table(cnee.res$phyloAcc_ra, cnee.res$phyloP_ra)
fisher.test(table(cnee.res$phyloAcc_ra, cnee.res$phyloP_ra))

#where to draw the line? look at posterior matrix and individual clade phyloP
ratite.tips<-c("aptRow", "aptOwe", "aptHaa", "casCas", "droNov", "rheAme", "rhePen", "anoDid", "strCam")
tips<- zpost %>% select(contains("_3")) %>% select(X0_3:X35_3)
tips[tips >= 0.90] = 1
tips[tips < 1] = 0
names(tips) = names(zmat)[1:36]
tips$total_losses = rowSums(tips)
tips$ratite_losses = rowSums(tips[ratite.tips])
tips$ID = lik$ID

cnee.res <- cnee.res %>% full_join(., tips, by=c("ID" = "ID"))
cnee.res$nonratite_losses = cnee.res$total_losses - cnee.res$ratite_losses

#ratite-accelerated CNEEs:
table(cnee.res$phyloAcc_ra) #18209

#now define these as ratite-specific and then as convergent
#defining ratite-specific, use Bayes factor
table(cnee.res$phyloAcc_ra, cnee.res$bf_full > 5)
table(cnee.res$phyloAcc_ra, cnee.res$nonratite_losses == 0)
table(cnee.res$phyloAcc_ra, cnee.res$nonratite_losses == 0 & cnee.res$bf_full > 5)

#ratite specific
cnee.res$phyloAcc_rspecific = ifelse(cnee.res$nonratite_losses == 0 & cnee.res$bf_full > 5, TRUE, FALSE)
table(cnee.res$phyloAcc_ra, cnee.res$phyloAcc_rspecific) #4189 ratite specific, ratite accelerated elements

#find convergences
convergences <- cnee.res %>% 
  mutate(rhea = as.logical(rheAme + rhePen), cas = as.logical(casCas + droNov), kiwi = as.logical(aptHaa + aptRow + aptOwe), ostrich = as.logical(strCam), moa = as.logical(anoDid), clade_ct = rhea + cas + kiwi + ostrich + moa) %>% 
  select(ID, rhea:clade_ct) %>% 
  full_join(., select(anoDid, ID=name, anoDid.phast = qval), by=c("ID" = "ID")) %>% 
  full_join(., select(strCam, ID=name, strCam.phast = qval), by=c("ID"="ID")) %>% 
  full_join(., select(kiwi, ID=name, kiwi.phast = qval), by=c("ID"="ID")) %>% 
  full_join(., select(rhea, ID=name, rhea.phast = qval), by=c("ID"="ID")) %>% 
  full_join(., select(casuar, ID=name, cas.phast = qval), by=c("ID"="ID")) %>% 
  mutate(clade_ct_phast = (anoDid.phast < 0.05) + (strCam.phast < 0.05) + (kiwi.phast < 0.05) + (rhea.phast < 0.05) + (cas.phast < 0.05))

cnee.res <- full_join(cnee.res, convergences, by=c("ID" = "ID"))

#clade counts for ratite-specific
table(cnee.res$clade_ct, cnee.res$phyloAcc_ra & cnee.res$phyloAcc_rspecific)
table(cnee.res$clade_ct_phast, cnee.res$phyloAcc_ra & cnee.res$phyloAcc_rspecific)


