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

#PLOTTING ORIGINAL##
gwide <- ext_gm %>% filter(rar_qval < 0.01, window_size == "1000kb_100kb_slide") %>% summarize(maxpval=max(rar_pval)) %>% pull(maxpval)
ss <- ext_gm %>% filter(rar_qval < 0.05, window_size == "1000kb_100kb_slide") %>% summarize(maxpval = max(rar_pval)) %>% pull(maxpval)

#convert to format for qqman
se_manh <- ext_gm %>% filter(window_size == "1000kb_100kb_slide") %>% select(chr, start, end, rar_pval, window) %>%
  mutate(chr = sub("chr", "", chr)) %>% 
  mutate(CHR = as.numeric(sub("Z", 29, chr))) %>%
  mutate(BP = start+((end+1-start)/2)) %>% 
  dplyr::rename(SNP = window) %>%
  dplyr::rename(P = rar_pval) %>% select(SNP, CHR, BP, P)


#write plot to pdf

pdf("fig3c.pdf")

manhattan(se_manh, chrlabs=c(1:15, 17:28, "Z"), suggestiveline=-log10(ss), genomewideline = -log10(gwide), col=c("blue4", "orange3"), type="l") 

dev.off()

#windows with an excess of convergent rars - manually annotate in figure
ext_gm %>% filter(window_size == "1000kb_100kb_slide") %>% mutate(RAR = rar_pval < gwide, CRAR = crar_pval < gwide, logp = -log10(rar_pval)) %>% filter(CRAR | RAR) %>% select(chr, start, end, RAR, CRAR, logp) %>% write_tsv("gwide_sig_windows.tsv")

