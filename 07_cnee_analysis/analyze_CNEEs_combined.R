library(tidyverse)

#master script to run phyloAcc and phyloP parsing and merge


#get phylop and phyloacc
source("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/analyze_phyloP.R")
source("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/analyze_phyloAcc.R")

#merge
cnee <- phyloacc %>% mutate(phylop.ratite = cnee %in% allRatite$name[allRatite$qval <= 0.05])
cnee <- cnee %>% mutate(phylop.kiwi = cnee %in% kiwi$name[kiwi$qval <= 0.05]) %>%
  mutate(phylop.moa = cnee %in% anoDid$name[anoDid$qval <= 0.05]) %>%
  mutate(phylop.cas = cnee %in% casuar$name[casuar$qval <= 0.05]) %>%
  mutate(phylop.ost = cnee %in% strCam$name[strCam$qval <= 0.05]) %>%
  mutate(phylop.rhea = cnee %in% rhea$name[rhea$qval <= 0.05]) %>%
  mutate(phylop.tin = cnee %in% tinamou$name[tinamou$qval <= 0.05]) %>%
  mutate(phylop.rcount = as.numeric(phylop.cas) + as.numeric(phylop.kiwi) + as.numeric(phylop.rhea) + as.numeric(phylop.moa) + as.numeric(phylop.ost))
  
#reminder: 
#ratite_loss = allow as many losses as are estimated by model
#ratite_loss_cons = max 1 loss per kiwi / rhea / emu cas / ostrich / moa 
#ratite_loss_cons_min = max 1 loss per kiwk rhea emu cas / ostrich / moa
# .prob suffix = uses raw posterior probabilities
# .mat suffix rounds < 0.95 to 0 and >= 0.95 to 1

#define properties
cnee <- cnee %>% mutate(ratite_accel.1 = bf1 >= 10 & ratite_loss.prob >= 1, ratite_accel.2 = ratite_accel.1 & phylop.ratite) %>%
  mutate(ratite_spec.1 = bf2 > 1 & nonratite_loss.prob < 1 & nonratite_loss.mat == 0, ratite_spec.2 = ratite_spec.1 & !phylop.tin) %>%
  mutate(ratite_conv.1 = ratite_loss_cons.prob > 2, ratite_conv.2 = ratite_conv.1 & ratite_loss_cons_min.mat >= 2 & phylop.rcount > 1)

#load annotation
gene<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/annotation/cnee/cnees.galgal4.annotation", header=F, stringsAsFactors = F) %>% tbl_df %>% rename(gene=V2)

#load bed
bed<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/annotation/cnee/cnees.galGal4UCSC.bed", stringsAsFactors = F) %>% tbl_df %>% rename(chr=V1, start=V2, end=V3, ID=V4)

cnee <- full_join(cnee, gene, by=c("cnee" = "V1")) %>% full_join(bed, by=c("cnee" = "ID"))

write_tsv(cnee, "~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/cnees.tsv") 
