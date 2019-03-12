## CODE TO ANALYZE CNEES ##
## UPDATED MARCH 2019 FOR FINAL REVISIONS AND TO ACCOMPANY DRYAD PACKAGE ##

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

#old cnees
final_v1 <- read_tsv("cnees.tsv.gz")

#define groups

define_groups <- function(x) {
  x %>%   
    mutate(accelerated = ifelse(logBF1 >= 10 & logBF2 > 1 & (it_pp_loss + ti_pp_loss) < 1, TRUE, FALSE)) %>%
    mutate(conv = ifelse(floss_cl_pp >= 1.80, TRUE, FALSE)) %>%
    mutate(phylop = as.logical(phylop)) %>% 
    mutate(dollo = ifelse(floss_cl_pp_dollo >= 1.80, TRUE, FALSE)) %>%
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
    mutate(conv = ifelse(floss_cl_pp >= 1.80, TRUE, FALSE)) %>%
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

#verions
final_original %>% define_groups() %>% group_by(version) %>% count(score) %>% spread(score,n)
final_reduced %>% define_groups() %>% group_by(version) %>% count(score) %>% spread(score,n)
final_extended %>% define_groups_corm() %>% group_by(version) %>% count(score) %>% spread(score,n)


##MAKE BEDS##

#original

pos_gg4 %>% 
  inner_join(final_original) %>% filter(version=="gain") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gapmissing_original_ucsc_galgal4.bed", col_names = FALSE)
  
pos_gg4 %>% 
  inner_join(final_original) %>% filter(version=="gain_gap") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gaploss_original_ucsc_galgal4.bed", col_names = FALSE)

pos_gg4_ncbi %>% 
  inner_join(final_original) %>% filter(version=="gain") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gapmissing_original_ncbi_galgal4.bed", col_names = FALSE)

pos_gg4_ncbi %>% 
  inner_join(final_original) %>% filter(version=="gain_gap") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gaploss_original_ncbi_galgal4.bed", col_names = FALSE)

pos_gg5 %>% 
  inner_join(final_original) %>% filter(version=="gain") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gapmissing_original_ucsc_galgal5.bed", col_names = FALSE)

pos_gg5 %>% 
  inner_join(final_original) %>% filter(version=="gain_gap") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gaploss_original_ucsc_galgal5.bed", col_names = FALSE)

pos_gg5_ncbi %>% 
  inner_join(final_original) %>% filter(version=="gain") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gapmissing_original_ncbi_galgal5.bed", col_names = FALSE)

pos_gg5_ncbi %>% 
  inner_join(final_original) %>% filter(version=="gain_gap") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gaploss_original_ncbi_galgal5.bed", col_names = FALSE)

#reduced

pos_gg4 %>% 
  inner_join(final_reduced) %>% filter(version=="gain") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gapmissing_reduced_ucsc_galgal4.bed", col_names = FALSE)

pos_gg4 %>% 
  inner_join(final_reduced) %>% filter(version=="gain_gap") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gaploss_reduced_ucsc_galgal4.bed", col_names = FALSE)

pos_gg4_ncbi %>% 
  inner_join(final_reduced) %>% filter(version=="gain") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gapmissing_reduced_ncbi_galgal4.bed", col_names = FALSE)

pos_gg4_ncbi %>% 
  inner_join(final_reduced) %>% filter(version=="gain_gap") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gaploss_reduced_ncbi_galgal4.bed", col_names = FALSE)

pos_gg5 %>% 
  inner_join(final_reduced) %>% filter(version=="gain") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gapmissing_reduced_ucsc_galgal5.bed", col_names = FALSE)

pos_gg5 %>% 
  inner_join(final_reduced) %>% filter(version=="gain_gap") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gaploss_reduced_ucsc_galgal5.bed", col_names = FALSE)

pos_gg5_ncbi %>% 
  inner_join(final_reduced) %>% filter(version=="gain") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gapmissing_reduced_ncbi_galgal5.bed", col_names = FALSE)

pos_gg5_ncbi %>% 
  inner_join(final_reduced) %>% filter(version=="gain_gap") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gaploss_reduced_ncbi_galgal5.bed", col_names = FALSE)

## extended

pos_gg4 %>% 
  inner_join(final_extended) %>% filter(version=="gain") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gapmissing_extended_ucsc_galgal4.bed", col_names = FALSE)

pos_gg4 %>% 
  inner_join(final_extended) %>% filter(version=="gain_gap") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gaploss_extended_ucsc_galgal4.bed", col_names = FALSE)

pos_gg4_ncbi %>% 
  inner_join(final_extended) %>% filter(version=="gain") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gapmissing_extended_ncbi_galgal4.bed", col_names = FALSE)

pos_gg4_ncbi %>% 
  inner_join(final_extended) %>% filter(version=="gain_gap") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gaploss_extended_ncbi_galgal4.bed", col_names = FALSE)

pos_gg5 %>% 
  inner_join(final_extended) %>% filter(version=="gain") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gapmissing_extended_ucsc_galgal5.bed", col_names = FALSE)

pos_gg5 %>% 
  inner_join(final_extended) %>% filter(version=="gain_gap") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gaploss_extended_ucsc_galgal5.bed", col_names = FALSE)

pos_gg5_ncbi %>% 
  inner_join(final_extended) %>% filter(version=="gain") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gapmissing_extended_ncbi_galgal5.bed", col_names = FALSE)

pos_gg5_ncbi %>% 
  inner_join(final_extended) %>% filter(version=="gain_gap") %>%
  define_groups %>% select(chr, start, end, cnee, score) %>%
  write_tsv(path="cnee_gaploss_extended_ncbi_galgal5.bed", col_names = FALSE)

## WRITE OUT FINAL DATA FILES

final_extended %>% define_groups_corm() %>% write_tsv("final_extended_cnee.tsv")
final_original %>% define_groups() %>% write_tsv("final_original_cnee.tsv")
final_reduced %>% define_groups() %>% write_tsv("final_reduced_cnee.tsv")


#table s8
final_original %>% filter(logBF1 >= 10, logBF2 >= 1, (ti_pp_loss + it_pp_loss) < 1, version=="gain") %>%
  select(cd_pp_loss, rh_pp_loss, os_pp_loss, ki_pp_loss, mo_pp_loss) %>% 
  mutate(cas_emu = cd_pp_loss > 0.80, rhea = rh_pp_loss > 0.80, ostrich = os_pp_loss > 0.80, kiwi = ki_pp_loss > 0.80, moa = mo_pp_loss > 0.80) %>% rowwise() %>%
  mutate(key = paste(cas_emu, rhea, kiwi, ostrich, moa, sep="-")) %>% count(key) %>% View()

final_reduced %>% filter(logBF1 >= 10, logBF2 >= 1, (ti_pp_loss + it_pp_loss) < 1, version=="gain") %>%
  select(cd_pp_loss, rh_pp_loss, os_pp_loss, ki_pp_loss) %>% 
  mutate(cas_emu = cd_pp_loss > 0.80, rhea = rh_pp_loss > 0.80, ostrich = os_pp_loss > 0.80, kiwi = ki_pp_loss > 0.80, moa = NA) %>% rowwise() %>%
  mutate(key = paste(cas_emu, rhea, kiwi, ostrich, moa, sep="-")) %>% count(key) %>% View()

