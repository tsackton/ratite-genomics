library(dplyr)

#load atac
atac <- read.table("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/atac/tissue_specific_R_May31.txt", header=T, stringsAsFactors = F) %>% tbl_df

#get phylop and phyloacc
source("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/analyze_phyloP.R")
source("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/analyze_phyloAcc.R")

#merge and clean up
atac <- atac %>% select(CNEE:Strn10_pool, St16_whole_H3K27ac:TotalRepressed)
cnee <- full_join(phyloacc, atac, by=c("cnee" = "CNEE"))
cnee <- cnee %>% mutate(phylop.ratite = cnee %in% allRatite$name[allRatite$qval <= 0.05])
cnee <- cnee %>% mutate(phylop.kiwi = cnee %in% kiwi$name[kiwi$qval <= 0.05]) %>%
  mutate(phylop.moa = cnee %in% anoDid$name[anoDid$qval <= 0.05]) %>%
  mutate(phylop.cas = cnee %in% casuar$name[casuar$qval <= 0.05]) %>%
  mutate(phylop.ost = cnee %in% strCam$name[strCam$qval <= 0.05]) %>%
  mutate(phylop.rhea = cnee %in% rhea$name[rhea$qval <= 0.05]) %>%
  mutate(phylop.tin = cnee %in% tinamou$name[tinamou$qval <= 0.05]) %>%
  mutate(phylop.rcount = as.numeric(phylop.cas) + as.numeric(phylop.kiwi) + as.numeric(phylop.rhea) + as.numeric(phylop.moa) + as.numeric(phylop.ost))
  

#define properties

cnee <- cnee %>% mutate(ratite_accel.1 = bf1 >= 10 & ratite_loss.prob >= 1, ratite_accel.2 = ratite_accel.1 & phylop.ratite) %>%
  mutate(ratite_spec.1 = bf2 > 1 & nonratite_loss.prob < 1 & nonratite_loss.mat == 0, ratite_spec.2 = ratite_spec.1 & !phylop.tin) %>%
  mutate(ratite_conv.1 = ratite_loss_cons.prob > 2, ratite_conv.2 = ratite_conv.1 & ratite_loss_cons_min.mat >= 2 & phylop.rcount > 1)

#load annotation
gene<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/annotation/cnee/cnees.galgal4.annotation", header=F, stringsAsFactors = F) %>% tbl_df %>% rename(gene=V2)

cnee <- full_join(cnee, gene, by=c("cnee" = "V1"))

##FILTERS FOR CANDIDATES###

#FL specific putative enhancers -- broad
cnee %>% filter(TotalStrict == 1, TotalActive >= 1, FL_strict > 0, ratite_accel.1, ratite_conv.1, ratite_spec.1) %>% select(cnee, FL_strict, TotalActive, TotalRepressed, ratite_loss_cons.mat, ratite_loss_cons.prob, bf1, bf2, gene, phylop.ratite, phylop.rcount) %>% arrange(ratite_loss_cons.mat, FL_strict) %>% print.data.frame

#FL specific putative enhancers -- narrow
cnee %>% filter(TotalStrict == 1, TotalActive >= 1, FL_strict > 0, ratite_accel.2, ratite_conv.1, ratite_spec.2) %>% select(cnee, FL_strict, TotalActive, TotalRepressed, ratite_loss_cons.mat, ratite_loss_cons.prob, bf1, bf2, gene,phylop.ratite, phylop.rcount) %>% arrange(ratite_loss_cons.mat, FL_strict) %>% print.data.frame

#FL specific putative enhancers -- good compromise
cnee %>% filter(TotalStrict == 1, TotalActive >= 1, FL_strict > 0, ratite_accel.2, ratite_spec.2, !phylop.tin, ratite_loss_cons_min.prob >= 1, ratite_loss_cons.mat >= 1) %>% select(cnee, gene, TotalActive, TotalRepressed, ratite_loss_cons.mat, ratite_loss_cons_min.prob, bf1, bf2, phylop.ratite, phylop.rcount) %>% arrange(gene) %>% print.data.frame

#Always on  putative enhancers -- good compromise
cnee %>% filter(TotalStrict >= 5, TotalActive >= 5, ratite_accel.2, ratite_spec.2, !phylop.tin, ratite_loss_cons_min.prob > 1.5, ratite_loss_cons.mat > 1) %>% select(cnee, gene, TotalStrict, TotalActive, TotalRepressed, ratite_loss_cons.mat, ratite_loss_cons_min.prob, bf1, bf2, phylop.ratite, phylop.rcount) %>% arrange(gene) %>% print.data.frame


#drop chip-seq - 13 elements
cnee %>% filter(AllBayes == 1, TotalStrict == 1, FL_strict > 0) %>% select(Chromosome:HL_strict, TotalActive, TotalRepressed) %>% arrange(Chromosome, Start) %>% print.data.frame

#relax convergent - 17 elements
cnee %>% filter(set4_conv1 > 0 | set4_conv2 > 0 | set4_conv3 > 0 | set4_conv4 > 0, TotalStrict == 1, EverActive > 0, FL_strict > 0, TotalRepressed+TotalActive <= TotalActive | TotalActive > 0) %>% select(Chromosome:HL_strict, TotalActive, TotalRepressed) %>% arrange(Chromosome, Start) %>% print.data.frame

#relax FL-specificity - 43 elements
cnee %>% filter(AllBayes == 1, TotalStrict > 0, EverActive > 0, FL_strict > 0, TotalRepressed+TotalActive <= TotalActive | TotalActive > 0) %>% select(Chromosome:HL_strict, TotalActive, TotalRepressed) %>% arrange(Chromosome, Start) %>% print.data.frame

#broadest set - 147 elements
cnee %>% filter(set4_conv1 > 0 | set4_conv2 > 0 | set4_conv3 > 0 | set4_conv4 > 0, FL_strict > 0, EverActive > 0, TotalRepressed+TotalActive <= TotalActive | TotalActive > 0) %>% select(Chromosome:HL_strict, TotalActive, TotalRepressed) %>% arrange(Chromosome, Start)
