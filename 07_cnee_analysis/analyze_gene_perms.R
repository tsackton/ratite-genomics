#revised Oct 2018 for paper revisions
library(tidyverse)

path_to_data <- "~/Projects/birds/ratite_compgen/DRYAD/07_cnees/processed"

#read real data
ext<-read_tsv(paste0(path_to_data, "/extended_gene_galgal4_run1_real.tsv.gz")) %>% rename(run = set) %>%
  mutate(total=in_target_TRUE + in_target_FALSE, 
         set = case_when(
           run == 1 ~ "rar",
           run == 2 ~ "crar",
           run == 3 ~ "crar_dollo")) %>%
  select(version, set, gene, count=in_target_TRUE, total)

orig<-read_tsv(paste0(path_to_data, "/original_gene_galgal4_run1_real.tsv.gz")) %>% rename(run = set) %>%
  mutate(total=in_target_TRUE + in_target_FALSE, 
         set = case_when(
           run == 1 ~ "rar",
           run == 2 ~ "crar",
           run == 3 ~ "crar_dollo")) %>%
  select(version, set, gene, count=in_target_TRUE, total)

red<-read_tsv(paste0(path_to_data, "/reduced_gene_galgal4_run1_real.tsv.gz")) %>% rename(run = set) %>%
  mutate(total=in_target_TRUE + in_target_FALSE, 
         set = case_when(
           run == 1 ~ "rar",
           run == 2 ~ "crar",
           run == 3 ~ "crar_dollo")) %>%
  select(version, set, gene, count=in_target_TRUE, total)


#read perms
ext_perm<-readRDS(paste0(path_to_data, "/geneperms/extended_galgal4.robj"))
orig_perm<-readRDS(paste0(path_to_data, "/geneperms/original_galgal4.robj"))
red_perm<-readRDS(paste0(path_to_data, "/geneperms/reduced_galgal4.robj"))

#compute P-value as 1-ecdf(real-1): ecdf(real-1) is prob (X <= (x-1)), e.g. prob(X < x), 1-that is prob (X >= x)

orig_merge <- full_join(orig, orig_perm, by=c("version" = "version", "set" = "set", "gene" = "gene")) %>% rowwise %>% mutate(pval = 1-ecdf_gene(count-1)) %>% ungroup %>% group_by(version, set) %>% mutate(pval = ifelse(pval == 0, 1e-04, pval))

ext_merge <- full_join(ext, ext_perm, by=c("version" = "version", "set" = "set", "gene" = "gene")) %>% rowwise %>% mutate(pval = 1-ecdf_gene(count-1)) %>% ungroup %>% group_by(version, set) %>% mutate(pval = ifelse(pval == 0, 1e-04, pval))

red_merge <- full_join(red, red_perm, by=c("version" = "version", "set" = "set", "gene" = "gene")) %>% rowwise %>% mutate(pval = 1-ecdf_gene(count-1)) %>% ungroup %>% group_by(version, set) %>% mutate(pval = ifelse(pval == 0, 1e-04, pval))

#plot
exclude_nonsig <- function(DF) {
  DF %>% filter(sig_class != "none")
}

spread_n <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

#clean up merged dataset for plotting
orig_merge_plot <- orig_merge %>%
  group_by(set) %>% 
  filter(version == "gain", count > 0) %>%
  mutate(qval = p.adjust(pval, "fdr")) %>%
  separate(gene, into=c("ncbi", "sym"), sep=":") %>%
  select(ncbi, sym, set, cnee_total = total, count, qval) %>%
  ungroup %>%
  spread_n(set, c(qval,count))

#change color scheme, change plotting character
orig_merge_plot %>% 
  mutate(sig_class = case_when(
    crar_qval <= 0.05 ~ "a2crar",
    rar_qval <= 0.05 ~ "a1rar",
    TRUE ~ "none"
  )) %>%
  ggplot(aes(x=rar_count, y=cnee_total, col=sig_class, label=sym)) + 
  theme_classic() + scale_y_log10() + geom_jitter(shape=16) +
  scale_color_brewer(palette="Dark2") +
  labs(x="Number of convergent RARs near gene", y="Total number of CNEEs near gene",color="Signficantly enriched?") +
  geom_text(data=exclude_nonsig, nudge_x = .1, show.legend=FALSE) + scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))

ggsave("Fig3B.pdf")

write_tsv(orig_merge_plot, "Fig3B-raw-data.tsv")
