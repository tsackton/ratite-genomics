## CODE TO ANALYZE CNEES ##
## UPDATED OCT 2018 FOR MANUSCRIPT REVISIONS ##

library(tidyverse)

#master script to run phyloAcc and phyloP parsing and merge

#get phylop and phyloacc
setwd("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/")
source("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/analyze_phyloAcc.R")
source("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/analyze_phyloP.R")
setwd("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/")

#add phlyoP results to each cnee set as appropriate

final_extended <- cnee_phyloP %>% select(name, contains("qval")) %>%
  mutate_at(vars(contains("qval")), ctb, cutoff=0.05, lower=TRUE) %>%
  rename_at(vars(contains("qval")), gsub, pattern="qval", replacement="phylop", fixed=TRUE) %>%
  full_join(cnee_ext_final, ., by=c("cnee" = "name")) %>% 
  mutate(floss_cl_phylop = phylop.mo + phylop.cd + phylop.ki + phylop.rh + phylop.os) %>%
  mutate(floss_cl_phylop_dollo = phylop.mo + tto(phylop.cd + phylop.ki + phylop.rh) + phylop.os)  

final_reduced <- cnee_phyloP %>% select(name, contains("qval")) %>%
  mutate_at(vars(contains("qval")), ctb, cutoff=0.05, lower=TRUE) %>%
  rename_at(vars(contains("qval")), gsub, pattern="qval", replacement="phylop", fixed=TRUE) %>%
  select(-phylop.mo) %>%
  full_join(cnee_red_final, ., by=c("cnee" = "name"))  %>% 
  mutate(floss_cl_phylop = phylop.cd + phylop.ki + phylop.rh + phylop.os) %>%
  mutate(floss_cl_phylop_dollo = tto(phylop.cd + phylop.ki + phylop.rh) + phylop.os)  

final_original <- cnee_phyloP %>% select(name, contains("qval")) %>%
  mutate_at(vars(contains("qval")), ctb, cutoff=0.05, lower=TRUE) %>%
  rename_at(vars(contains("qval")), gsub, pattern="qval", replacement="phylop", fixed=TRUE) %>%
  full_join(cnee_orig_final, ., by=c("cnee" = "name")) %>% 
  mutate(floss_cl_phylop = phylop.mo + phylop.cd + phylop.ki + phylop.rh + phylop.os) %>%
  mutate(floss_cl_phylop_dollo = phylop.mo + tto(phylop.cd + phylop.ki + phylop.rh) + phylop.os)

#load annotation
gene_gg4 <-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/04_wga/03_ce_annotation/cnees.galgal4.annotation", col_names=c("cnee", "gene"))
gene_gg5 <-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/04_wga/03_ce_annotation/cnees.galgal5.annotation", col_names=c("cnee", "gene"))
pos_gg4<-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/04_wga/03_ce_annotation/cnees.galGal4UCSC.bed", col_names = c("chr", "start", "end", "cnee"))
pos_gg4_ncbi<-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/04_wga/03_ce_annotation/cnees.galGal4NCBI.bed", col_names = c("chr", "start", "end", "cnee"))
pos_gg5<-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/04_wga/03_ce_annotation/cnees.galGal5UCSC.bed", col_names = c("chr", "start", "end", "cnee"))
pos_gg5_ncbi<-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/04_wga/03_ce_annotation/cnees.galGal5NCBI.bed", col_names = c("chr", "start", "end", "cnee"))

## old definitions of categories
#cnee <- cnee %>% mutate(ratite_accel.1 = bf1 >= 10 & ratite_loss.prob >= 1, ratite_accel.2 = ratite_accel.1 & phylop.ratite) %>%
#  mutate(ratite_spec.1 = bf2 > 1 & nonratite_loss.prob < 1 & nonratite_loss.mat == 0, ratite_spec.2 = ratite_spec.1 & !phylop.tin) %>%
#  mutate(ratite_conv.1 = ratite_loss_cons.prob > 2, ratite_conv.2 = ratite_conv.1 & ratite_loss_cons_min.mat >= 2 & phylop.rcount > 1)



#define groups

define_groups <- function(x) {
  x %>%   
    mutate(accelerated = ifelse(logBF1 >= 10 & logBF2 > 1 & (it_pp_loss + ti_pp_loss) < 1, TRUE, FALSE)) %>%
    mutate(conv = ifelse(floss_cl_pp >= 2, TRUE, FALSE)) %>%
    mutate(phylop = as.logical(phylop)) %>% 
    mutate(dollo = ifelse(floss_cl_pp_dollo >= 2, TRUE, FALSE)) %>%
    mutate(score = case_when(
      accelerated & conv & phylop & dollo ~ "Dollo_Conv_Acc_PhyloP",
      accelerated & conv & !phylop & !dollo ~ "Conv_Acc",
      accelerated & conv & phylop & !dollo ~ "Conv_Acc_PhyloP",
      accelerated & conv & !phylop & dollo ~ "Dollo_Conv_Acc",
      accelerated & !conv & phylop ~ "Acc_PhyloP",
      accelerated & !conv & !phylop ~ "Acc",
      TRUE ~ "none"
    ))
}

define_groups_corm <- function(x) {
  x %>%
    mutate(accelerated = ifelse(logBF1 >= 10 & logBF2 >= 1 & (ti_pp_loss + it_pp_loss) < 1, TRUE, FALSE)) %>%
    mutate(conv = ifelse(floss_cl_pp >= 2, TRUE, FALSE)) %>%
    mutate(ratite = ifelse((floss_cl_pp - gc_pp_loss) >= 0.90, TRUE, FALSE)) %>%
    mutate(corm = ifelse(gc_pp_loss >= 0.90, TRUE, FALSE)) %>%
    mutate(score = case_when(
      accelerated & conv & ratite & corm ~ "Conv_Acc_Both",
      accelerated & conv & ratite & !corm ~ "Conv_Acc_Ratite",
      accelerated & conv & !ratite & corm ~ "Conv_Acc_Corm",
      accelerated & conv & !ratite & !corm ~ "Conv_Acc_Neither",
      accelerated & !conv & ratite & corm ~ "Acc_Both",
      accelerated & !conv & ratite & !corm ~ "Acc_Ratite",
      accelerated & !conv & !ratite & corm ~ "Acc_Corm",
      accelerated & !conv & !ratite & !corm ~ "Acc_Neither",
      TRUE ~ "none"
    ))
}

pos_gg4 %>% 
  inner_join(final_original) %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_scored_original_ucsc_galgal4.bed", col_names = FALSE)
  
pos_gg4_ncbi %>% 
  inner_join(final_original) %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_scored_original_ncbi_galgal4.bed", col_names = FALSE)

pos_gg4 %>% 
  inner_join(final_reduced) %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_scored_reduced_ucsc_galgal4.bed", col_names = FALSE)

pos_gg4_ncbi %>% 
  inner_join(final_reduced) %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_scored_reduced_ncbi_galgal4.bed", col_names = FALSE)

pos_gg4 %>% 
  inner_join(final_extended) %>%
  define_groups_corm %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_scored_extended_ucsc_galgal4.bed", col_names = FALSE)

pos_gg4_ncbi %>% 
  inner_join(final_extended) %>%
  define_groups_corm %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_scored_extended_ncbi_galgal4.bed", col_names = FALSE)



pos_gg5 %>% 
  inner_join(final_original) %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_scored_original_ucsc_galgal5.bed", col_names = FALSE)

pos_gg5_ncbi %>% 
  inner_join(final_original) %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_scored_original_ncbi_galgal5.bed", col_names = FALSE)

pos_gg5 %>% 
  inner_join(final_reduced) %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_scored_reduced_ucsc_galgal5.bed", col_names = FALSE)

pos_gg5_ncbi %>% 
  inner_join(final_reduced) %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_scored_reduced_ncbi_galgal5.bed", col_names = FALSE)

pos_gg5 %>% 
  inner_join(final_extended) %>%
  define_groups_corm %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_scored_extended_ucsc_galgal5.bed", col_names = FALSE)

pos_gg5_ncbi %>% 
  inner_join(final_extended) %>%
  define_groups_corm %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_scored_extended_ncbi_galgal5.bed", col_names = FALSE)


#more testing
final_extended %>% define_groups %>% inner_join(gene_gg4) %>% filter(grepl("TBX5", gene)) %>% count(score)

#old cnees
final_v1 <- read_tsv("cnees.tsv.gz")

#compare
final_original %>% define_groups %>% full_join(final_v1, by=c("cnee" = "cnee")) %>%
  mutate(old_score = case_when(
    ratite_accel.1 & ratite_spec.1 & ratite_conv.1 ~ "Conv_Acc",
    ratite_accel.1 & ratite_spec.1 & !ratite_conv.1 ~ "Acc",
    TRUE ~ "none"
  )) %>% 
  mutate(score = sub("_PhyloP", "", score)) %>% 
  mutate(score = sub("Dollo_", "", score)) %>%
  with(., table(score, old_score))

final_original %>% define_groups %>% full_join(final_v1, by=c("cnee" = "cnee")) %>%
  mutate(old_score = case_when(
    ratite_accel.1 & ratite_spec.1 & ratite_conv.1 ~ "Conv_Acc",
    ratite_accel.1 & ratite_spec.1 & !ratite_conv.1 ~ "Acc",
    TRUE ~ "none"
  )) %>% 
  mutate(score = sub("_PhyloP", "", score)) %>% 
  mutate(score = sub("Dollo_", "", score)) %>%
  select(cnee:l_rate, bf1:loglik.full, cd_pp_loss:floss_cl_phylop_dollo, cas_loss.prob:neo_loss.prob, ratite_accel.1:gene, score, old_score) %>%
  filter(score == "none" & old_score=="Conv_Acc") %>% View()


final_original %>% define_groups %>% full_join(final_v1, by=c("cnee" = "cnee")) %>% ggplot(aes(logBF2, bf2)) + geom_point()
final_original %>% define_groups %>% mutate(score = sub("_PhyloP", "", score)) %>% 
  mutate(score = sub("Dollo_", "", score)) %>%
  full_join(final_v1, by=c("cnee" = "cnee")) %>% filter(pmax(bf1,logBF1) > 10, pmax(bf2,logBF2)>1) %>% 
  ggplot(aes(logBF1, bf1)) + geom_point(alpha=0.1) + geom_hline(yintercept=10, col="red") + geom_vline(xintercept=10, col="red")+ coord_cartesian(xlim=c(-25,50), ylim=c(-25,50))

final_original %>% define_groups %>% full_join(final_v1, by=c("cnee" = "cnee")) %>% mutate(score = sub("_PhyloP", "", score)) %>% 
  mutate(score = sub("Dollo_", "", score)) %>%
  filter(pmax(bf1,logBF1) > 10) %>% ggplot(aes(floss_cl_pp_dollo, ratite_loss_cons_min.prob)) + geom_point(col="red", alpha=0.2) + geom_abline() 

          
final_original %>% full_join(final_v1, by=c("cnee" = "cnee")) %>% filter(ratite_loss_cons.prob > 4, floss_cl_pp < 1) %>% select(cnee, gene)

final_original %>% full_join(final_reduced, by=c("cnee" = "cnee")) %>% ggplot(aes(strCam.x, strCam.y)) + geom_point(alpha=0.1) + geom_abline(col="red", size=3, linetype="dashed")

final_original %>% full_join(final_reduced, by=c("cnee" = "cnee")) %>% ggplot(aes(strCam.x, strCam.y)) + geom_point(alpha=0.1) + geom_abline(col="red", size=3, linetype="dashed")

final_original %>% full_join(final_extended, by=c("cnee" = "cnee")) %>% ggplot(aes(strCam.x, strCam.y)) + geom_point(alpha=0.1) + geom_abline(col="red", size=3, linetype="dashed")

final_extended %>% full_join(final_original, by=c("cnee" = "cnee")) %>% ggplot(aes(logBF1.x, logBF1.y)) + geom_point(alpha=0.1) + coord_cartesian(xlim=c(-25,50), ylim=c(-25,50))

final_original %>% full_join(final_v1, by=c("cnee" = "cnee")) %>% ggplot(aes(strCam, strCam.prob)) + geom_point(alpha=0.1) + geom_abline(col="red", size=1, linetype="dashed")

final_original %>% full_join(final_v1, by=c("cnee" = "cnee")) %>% ggplot(aes(anoDid, anoDid.prob)) + geom_point(alpha=0.1) + geom_abline(col="red", size=1, linetype="dashed")

final_original %>% define_groups %>% full_join(final_v1, by=c("cnee" = "cnee")) %>%
  mutate(old_score = case_when(
    ratite_accel.1 & ratite_spec.1 & ratite_conv.1 ~ "Acc",
    ratite_accel.1 & ratite_spec.1 & !ratite_conv.1 ~ "Acc",
    TRUE ~ "none"
  )) %>% 
  mutate(score = sub("_PhyloP", "", score)) %>% 
  mutate(score = sub("Dollo_", "", score)) %>%
  mutate(score = sub("Conv_", "", score)) %>%
  filter(!(score == old_score)) %>%
  mutate(loss_diff = round(ratite_loss_cons_min.prob - floss_cl_pp_dollo,2)) %>%
  mutate(old_new = paste0(old_score, "_", score)) %>% 
  ggplot(aes((it_pp_loss+ti_pp_loss), (internal_loss.prob + tin_loss.prob), color=old_new)) + geom_point(alpha=0.4) + geom_abline(col="red", size=1, linetype="dashed")

final_reduced %>% full_join(final_v1, by=c("cnee" = "cnee")) %>%
  mutate(bfAll_NEW = (loglik_Full - loglik_Null)) %>%
  mutate(bfAll_OLD = (loglik.full - loglik.null)) %>% 
  ggplot(aes(bfAll_NEW, bfAll_OLD)) + geom_point(alpha=0.1) + coord_cartesian(xlim=c(-30,75), ylim=c(-30,75))

