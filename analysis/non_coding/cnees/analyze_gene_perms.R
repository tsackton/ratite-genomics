
library(data.table)
library(tidyverse)
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


obs_gene_pval %>% filter(epval <= 1e-03, set=="set3") %>% select(gene, in_target_TRUE, in_target_FALSE) %>% print.data.frame
obs_gene_pval %>% filter(epval <= 1e-03, set=="set1") %>% select(gene) %>% print.data.frame

obs_gene_pval %>% filter(set=="set3") %>% mutate(cnee_total = in_target_FALSE + in_target_TRUE, sig = epval < 1e-03) %>% ggplot(aes(x=in_target_TRUE, y=cnee_total, col=sig)) + theme_classic() + scale_y_log10() +  geom_jitter()
obs_gene_pval %>% filter(set=="set1") %>% mutate(cnee_total = in_target_FALSE + in_target_TRUE, sig = epval < 1e-03) %>% ggplot(aes(x=in_target_TRUE, y=cnee_total, col=sig)) + scale_y_log10() +  geom_jitter()
