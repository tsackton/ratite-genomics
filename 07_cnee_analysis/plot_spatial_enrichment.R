#plot spatial enrichment

library("tidyverse")
library(qqman)

se_results <- read_tsv("spatial_enrichment.results", guess_max = 10000)

se_results %>% filter(rar.1.qval < 0.01, window_size == "1000kb_100kb_slide") %>% arrange(chr, start) %>% select(window_size, chr, start, end, rar.1, rar.1.pval) %>%print.data.frame
se_results %>% filter(crar.1.qval < 0.1, window_size == "1000kb_100kb_slide") %>% arrange(chr, start) %>% select(window_size, chr, start, end, crar.1, crar.1.pval) %>%print.data.frame

#pval for q = 0.05
se_results %>% filter(crar.1.qval < 0.01, window_size == "1000kb_100kb_slide") %>% summarize(max(crar.1.pval))
se_results %>% filter(crar.1.qval < 0.05, window_size == "1000kb_100kb_slide") %>% summarize(max(crar.1.pval))

#convert to format for qqman
se_manh <- se_results %>% filter(window_size == "1000kb_100kb_slide") %>% select(chr, start, end, crar.1.pval, window) %>%
  mutate(chr = sub("chr", "", chr)) %>% 
  mutate(CHR = as.numeric(sub("Z", 29, chr))) %>%
  mutate(BP = start+((end+1-start)/2)) %>% 
  dplyr::rename(SNP = window) %>%
  dplyr::rename(P = crar.1.pval) %>% select(SNP, CHR, BP, P)

pdf("~/Projects/birds/ratite_compgen/manuscript/ver7/Figure3E_2_2.pdf")
manhattan(se_manh, chrlabs=c(1:15, 17:24, 26, 27, 28, "Z"), suggestiveline=-log10(0.000552), genomewideline = -log10(0.0000654), type="l", col=c("blue4", "orange3")) 
dev.off()

se_results %>% filter(crar.1.qval < 0.01, window_size == "1000kb_100kb_slide") %>% select(chr, start, end) %>% print.data.frame

##NOT USED##
#genes
se_results %>% filter(chr=="chr6", window_size == "1000kb_100kb_slide", crar.1.qval < 0.1) %>% select(start, end) %>% print.data.frame


se_results %>% filter(window_size == "1000kb_100kb_slide") %>% 
  mutate(logp = -1 * log10(crar.1.pval)) %>% 
  ggplot(aes(x=window, y=logp)) + theme_classic() + geom_line() + 
  labs(x="Window", y="-log10(pval) for enrichment of convergent RARs") +
  geom_hline(yintercept = -log10(0.000552), col="red",linetype="dashed") + 
  annotate("text", y=-log10(0.000552)+0.2, x=500, label="5% FDR")
ggsave("~/Projects/birds/ratite_compgen/manuscript/ver7/Figure3D_2.pdf")


se_results %>% filter(chr=="chr7", window_size == "1000kb_100kb_slide") %>% 
  mutate(logp = -1 * log10(crar.1.pval)) %>% 
  ggplot(aes(x=start, y=logp)) + theme_classic() + geom_line() + 
  labs(x="Position on chromosome 7", y="-log10(pval) for enrichment of convergent RARs") +
  geom_hline(yintercept = -log10(0.000552), col="red",linetype="dashed") + 
  annotate("text", y=-log10(0.000552)+0.2, x=0, label="5% FDR")
ggsave("~/Projects/birds/ratite_compgen/manuscript/ver5/Figure3E.pdf")
