#plot spatial enrichment - clean up to do right

library(tidyverse)
library(qqman)

pdf("~/Downloads/spatial_plots_2.pdf")
se_results <- read_tsv("extended_spatial_results.tsv", guess_max = 10000) %>% filter(version=="gap")
se_results <- read_tsv("extended_spatial_results.tsv", guess_max = 10000) %>% filter(version=="gain_gap")
se_results <- read_tsv("extended_spatial_results.tsv", guess_max = 10000) %>% filter(version=="orig")
se_results <- read_tsv("extended_spatial_results.tsv", guess_max = 10000) %>% filter(version=="gain")

#pval for q = 0.05
gwide <- se_results %>% filter(crar_qval < 0.01, window_size == "1000kb_100kb_slide") %>% summarize(maxpval=max(crar_pval)) %>% pull(maxpval)
ss <- se_results %>% filter(crar_qval < 0.05, window_size == "1000kb_100kb_slide") %>% summarize(maxpval = max(crar_pval)) %>% pull(maxpval)

#convert to format for qqman
se_manh <- se_results %>% filter(window_size == "1000kb_100kb_slide") %>% select(chr, start, end, crar_pval, window) %>%
  mutate(chr = sub("chr", "", chr)) %>% 
  mutate(CHR = as.numeric(sub("Z", 29, chr))) %>%
  mutate(BP = start+((end+1-start)/2)) %>% 
  dplyr::rename(SNP = window) %>%
  dplyr::rename(P = crar_pval) %>% select(SNP, CHR, BP, P)


manhattan(se_manh, chrlabs=c(1:15, 17:24, 26, 27, 28, "Z"), suggestiveline=-log10(0.000204), genomewideline = -log10(0.0000278), type="l", col=c("blue4", "orange3")) 

dev.off()