#plot spatial enrichment - clean up to do right

library(tidyverse)
library(qqman)
setwd("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/")
orig_gm <- read_tsv("original_spatial_results.tsv.gz", guess_max = 10000) %>% filter(version=="gain")
orig_gl <- read_tsv("original_spatial_results.tsv.gz", guess_max = 10000) %>% filter(version=="gain_gap")
ext_gm <- read_tsv("extended_spatial_results.tsv.gz", guess_max = 10000) %>% filter(version=="gain")
ext_gl <- read_tsv("extended_spatial_results.tsv.gz", guess_max = 10000) %>% filter(version=="gain_gap")
red_gm <- read_tsv("reduced_spatial_results.tsv.gz", guess_max = 100000) %>% filter(version=="gain")
red_gl <- read_tsv("reduced_spatial_results.tsv.gz", guess_max = 100000) %>% filter(version=="gain_gap")
corm_gm <- read_tsv("extended_ratiteVcorm_spatial_results.tsv.gz", guess_max = 10000) %>% filter(version=="gain")
corm_gl <- read_tsv("extended_ratiteVcorm_spatial_results.tsv.gz", guess_max = 10000) %>% filter(version=="gain_gap")

#PLOTTING ORIGINAL GAP##
gwide <- ext_gm %>% filter(rar_qval < 0.01, window_size == "1000kb_100kb_slide") %>% summarize(maxpval=max(rar_pval)) %>% pull(maxpval)
ss <- ext_gm %>% filter(rar_qval < 0.05, window_size == "1000kb_100kb_slide") %>% summarize(maxpval = max(rar_pval)) %>% pull(maxpval)

#convert to format for qqman
se_manh <- ext_gm %>% filter(window_size == "1000kb_100kb_slide") %>% select(chr, start, end, rar_pval, window) %>%
  mutate(chr = sub("chr", "", chr)) %>% 
  mutate(CHR = as.numeric(sub("Z", 29, chr))) %>%
  mutate(BP = start+((end+1-start)/2)) %>% 
  dplyr::rename(SNP = window) %>%
  dplyr::rename(P = rar_pval) %>% select(SNP, CHR, BP, P)


manhattan(se_manh, chrlabs=c(1:15, 17:28, "Z"), suggestiveline=-log10(ss), genomewideline = -log10(gwide), col=c("blue4", "orange3"), type="l") 

#windows with an excess of convergent rars - manually annotate in figure
ext_gm %>% filter(window_size == "1000kb_100kb_slide") %>% mutate(RAR = rar_pval < ss, CRAR = crar_pval < ss, logp = -log10(rar_pval)) %>% filter(CRAR) %>% select(chr, start, end, RAR, CRAR, logp) %>% View()

#location of cormorant accelerated elements
corm_elem<-read_tsv("final_extended_cnee.tsv") %>% filter(gc_pp_loss > 0.9, logBF1 >= 10, logBF2 >= 1)
cnee_pos <- read_tsv("cnee_gaploss_reduced_ucsc_galgal4.bed", col_names = c("chr", "start", "end", "cnee", "score")) %>% 
  mutate(pos = round((start+end)/2,0)) %>%
  select(chr, cnee, pos) 
corm_elem %>% filter(version == "gain") %>% select(cnee, gc_pp_loss, cd_pp_loss, rh_pp_loss, os_pp_loss, ki_pp_loss, mo_pp_loss) %>% left_join(cnee_pos, by=c("cnee" = "cnee")) %>% arrange(chr, pos) %>% print.data.frame

