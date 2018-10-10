## CODE TO ANALYZE CNEES ##
## UPDATED OCT 2018 FOR MANUSCRIPT REVISIONS ##

library(tidyverse)

#master script to run phyloAcc and phyloP parsing and merge

#get phylop and phyloacc
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

final_ver1 <- cnee_phyloP %>% select(name, contains("qval")) %>%
  mutate_at(vars(contains("qval")), ctb, cutoff=0.05, lower=TRUE) %>%
  rename_at(vars(contains("qval")), gsub, pattern="qval", replacement="phylop", fixed=TRUE) %>%
  full_join(cnee_v1_final, ., by=c("cnee" = "name")) %>% 
  mutate(floss_cl_phylop = phylop.mo + phylop.cd + phylop.ki + phylop.rh + phylop.os) %>%
  mutate(floss_cl_phylop_dollo = phylop.mo + tto(phylop.cd + phylop.ki + phylop.rh) + phylop.os)


#load annotation
gene_gg4 <-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/04_wga/03_ce_annotation/cnees.galgal4.annotation", col_names=c("cnee", "gene"))
gene_gg5 <-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/04_wga/03_ce_annotation/cnees.galgal5.annotation", col_names=c("cnee", "gene"))
pos_gg4<-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/04_wga/03_ce_annotation/cnees.galGal4UCSC.bed", col_names = c("chr", "start", "end", "cnee"))
pos_gg4_ncbi<-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/04_wga/03_ce_annotation/cnees.galGal4NCBI.bed", col_names = c("chr", "start", "end", "cnee"))
pos_gg5<-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/04_wga/03_ce_annotation/cnees.galGal5UCSC.bed", col_names = c("chr", "start", "end", "cnee"))
pos_gg5_ncbi<-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/04_wga/03_ce_annotation/cnees.galGal5NCBI.bed", col_names = c("chr", "start", "end", "cnee"))


#define groups

define_groups <- function(x) {
  x %>%   
    mutate(accelerated = ifelse(bf1 >= 10 & bf2 >= 1, TRUE, FALSE)) %>%
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
    mutate(accelerated = ifelse(bf1 >= 10 & bf2 >= 1, TRUE, FALSE)) %>%
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



final_ver1 %>% mutate(accelerated = ifelse(bf1 > 10, "Acc", "Cons")) %>% 
  mutate(nonratite_loss = ti_pp_loss + it_pp_loss) %>% filter(accelerated == "Acc") %>%
  mutate(ratite_specific = ifelse(bf2 > 1, "yes", "no")) %>%
  ggplot(aes(nonratite_loss, bf2, col=ratite_specific)) + 
  scale_color_brewer(palette="Set1") + 
  geom_point(alpha=0.2) + 
  guides(color=guide_legend(override.aes = list(alpha=1)))

final_ver1 %>% mutate(rar = ifelse(bf1 > 10, "Acc", "Cons")) %>%
  mutate(nonratite_loss = ifelse((nv_pp_loss + ti_pp_loss + it_pp_loss) > 1, "yes", "no")) %>%
  ggplot(aes(x=rar, y=bf2, fill=nonratite_loss)) + geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits=c(-20,20))
  
final_ver1 %>% mutate(rar = ifelse(bf1 > 10 & bf2 > 1, "Acc", "Cons")) %>%
  mutate(nonratite_loss = ifelse((nv_pp_loss + ti_pp_loss + it_pp_loss) > 1, "yes", "no")) %>%
  group_by(rar) %>%
  count(nonratite_loss)
