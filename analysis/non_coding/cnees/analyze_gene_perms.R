
library(data.table)
library(tidyverse)
perms<-fread("perm_gene_count_results.tsv")
setkey(perms, set, gene)
real_data<-fread("obs_gene_count_results.tsv") %>% mutate(set = paste0("set", set)) %>% as.tibble

get_empirical_pval <- function(gene_test, set_test, value) {
  null <- perms[.(set_test, gene_test)] %>% pull(rand_TRUE)
  (sum(null >= value)+1) / (length(null)+1)
}

obs_gene_pval <- real_data %>% group_by(set, gene) %>% mutate(epval = get_empirical_pval(gene, set, in_target_TRUE))
write_tsv(obs_gene_pval, "obs_gene_pval.results")