#revised Oct 2018 for paper revisions
library(data.table)
library(tidyverse)

setwd("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/")

#testing

file<-"geneperms/original_gene_galgal4_run69_perm.tsv"
perms<-read_tsv(file)
perms_ecdf <- perms %>% group_by(version, set, gene) %>% summarize(ecdf_gene = list(ecdf(rand_TRUE)))


perms_ecdf %>% filter(version == "gain", set == "rar", gene=="100216000:CNN2") %>% mutate(count = 5, pval=1-ecdf_gene[[1]](count))
#read permutation results

#this should work. then compute P-value as 1-ecdf(real-1)

perms<-fread("perm_gene_count_results.tsv")
setkey(perms, set, gene)
real_data<-fread("obs_gene_count_results.tsv") %>% mutate(set = paste0("set", set)) %>% as.tibble

get_empirical_pval <- function(gene_test, set_test, value) {
  null <- perms[.(set_test, gene_test)] %>% pull(rand_TRUE)
  (sum(null >= value)+1) / (length(null)+1)
}

#run
obs_gene_pval <- real_data %>% group_by(set, gene) %>% mutate(epval = get_empirical_pval(gene, set, in_target_TRUE))
write_tsv(obs_gene_pval, "obs_gene_pval.results")

#read
obs_gene_pval <- read_tsv("obs_gene_pval.results")


obs_gene_pval %>% filter(epval <= 1e-03, set=="set3") %>% dplyr::select(gene, in_target_TRUE, in_target_FALSE) %>% print.data.frame
obs_gene_pval %>% filter(epval <= 1e-03, set=="set1") %>% dplyr::select(gene) %>% print.data.frame

exclude_nonsig <- function(DF) {
  DF %>% filter(sig == TRUE)
}

obs_gene_pval %>% filter(set=="set3") %>% filter(in_target_TRUE > 0) %>%
  mutate(qval = p.adjust(epval, "fdr")) %>%
  mutate(cnee_total = in_target_FALSE + in_target_TRUE, sig = qval < 0.05) %>% 
  separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
  ggplot(aes(x=in_target_TRUE, y=cnee_total, col=sig, label=sym)) + theme_classic() + scale_y_log10() + geom_jitter() +
  labs(x="Number of convergent RARs near gene", y="Total number of CNEEs near gene",color="Signficantly enriched?") +
  geom_text(data=exclude_nonsig, check_overlap=FALSE, nudge_x = .3, show.legend=FALSE)
ggsave("~/Projects/birds/ratite_compgen/manuscript/ver7/Figure3B_2.pdf")

obs_gene_pval %>% filter(set=="set1")  %>% filter(in_target_TRUE > 0) %>%
  mutate(qval = p.adjust(epval, "fdr")) %>%
  mutate(cnee_total = in_target_FALSE + in_target_TRUE, sig = qval < 0.05) %>% 
  separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
  ggplot(aes(x=in_target_TRUE, y=cnee_total, col=sig, label=sym)) + theme_classic() + scale_y_log10() + geom_jitter() +
  labs(x="Number of RARs near gene", y="Total number of CNEEs near gene",color="Signficantly enriched?") +
  geom_text(data=exclude_nonsig, check_overlap=FALSE, nudge_x = .3, show.legend=FALSE)
ggsave("~/Projects/birds/ratite_compgen/manuscript/ver7/Figure3C_2.pdf")

