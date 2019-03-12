#revised Oct 2018 for paper revisions
library(tidyverse)

#read real data
ext<-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/geneperms/extended_gene_galgal4_run1_real.tsv") %>% rename(run = set) %>%
  mutate(total=in_target_TRUE + in_target_FALSE, 
         set = case_when(
           run == 1 ~ "rar",
           run == 2 ~ "crar",
           run == 3 ~ "crar_dollo")) %>%
  select(version, set, gene, count=in_target_TRUE, total)

orig<-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/geneperms/original_gene_galgal4_run1_real.tsv") %>% rename(run = set) %>%
  mutate(total=in_target_TRUE + in_target_FALSE, 
         set = case_when(
           run == 1 ~ "rar",
           run == 2 ~ "crar",
           run == 3 ~ "crar_dollo")) %>%
  select(version, set, gene, count=in_target_TRUE, total)

red<-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/geneperms/reduced_gene_galgal4_run1_real.tsv") %>% rename(run = set) %>%
  mutate(total=in_target_TRUE + in_target_FALSE, 
         set = case_when(
           run == 1 ~ "rar",
           run == 2 ~ "crar",
           run == 3 ~ "crar_dollo")) %>%
  select(version, set, gene, count=in_target_TRUE, total)


#read perms
ext_perm<-readRDS("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/geneperms/extended_galgal4.robj")
orig_perm<-readRDS("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/geneperms/original_galgal4.robj")
red_perm<-readRDS("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/geneperms/reduced_galgal4.robj")


#compute P-value as 1-ecdf(real-1): ecdf(real-1) is prob (X <= (x-1)), e.g. prob(X < x), 1-that is prob (X >= x)

orig_merge <- full_join(orig, orig_perm, by=c("version" = "version", "set" = "set", "gene" = "gene")) %>% rowwise %>% mutate(pval = 1-ecdf_gene(count-1)) %>% ungroup %>% group_by(version, set) %>% mutate(pval = ifelse(pval == 0, 1e-04, pval))

ext_merge <- full_join(ext, ext_perm, by=c("version" = "version", "set" = "set", "gene" = "gene")) %>% rowwise %>% mutate(pval = 1-ecdf_gene(count-1)) %>% ungroup %>% group_by(version, set) %>% mutate(pval = ifelse(pval == 0, 1e-04, pval))

red_merge <- full_join(red, red_perm, by=c("version" = "version", "set" = "set", "gene" = "gene")) %>% rowwise %>% mutate(pval = 1-ecdf_gene(count-1)) %>% ungroup %>% group_by(version, set) %>% mutate(pval = ifelse(pval == 0, 1e-04, pval))

#plot
exclude_nonsig <- function(DF, cutoff = 0.05) {
  DF %>% filter(qval <= cutoff)
}

ext_merge %>% filter(set=="rar", version=="gain") %>% filter(count > 0) %>%
  mutate(qval = p.adjust(pval, "fdr")) %>%
  mutate(cnee_total = total, sig = qval < 0.05) %>% 
  separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
  ggplot(aes(x=count, y=cnee_total, col=sig, label=sym)) + theme_classic() + scale_y_log10() + geom_jitter() +
  labs(x="Number of convergent RARs near gene", y="Total number of CNEEs near gene",color="Signficantly enriched?") +
  geom_text(data=exclude_nonsig, nudge_x = .1, show.legend=FALSE) + scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))

#also convergent RARs
orig_merge %>% filter(set=="rar", version=="gain") %>% filter(count > 0) %>%
  mutate(qval = p.adjust(pval, "fdr")) %>%
  mutate(cnee_total = total, sig = qval < 0.05) %>% 
  separate(gene, into=c("ncbi", "sym"), sep=":") %>% filter(qval < 0.05) %>%
  select(sym,count,total) %>% print.data.frame

#also convergent RARs
ext_merge %>% filter(set=="rar", version=="gain") %>% filter(count > 0) %>%
  mutate(qval = p.adjust(pval, "fdr")) %>%
  mutate(cnee_total = total, sig = qval < 0.05) %>% 
  separate(gene, into=c("ncbi", "sym"), sep=":") %>% filter(qval < 0.05) %>%
  select(sym,count,total) %>% print.data.frame
