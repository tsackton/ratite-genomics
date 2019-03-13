#revised Nov 2018 for paper revisions
library(tidyverse)

path_to_data <- "~/Projects/birds/ratite_compgen/DRYAD/07_cnees/processed"

transform_real <- function(DF) {
  DF %>% 
    separate(GeneRatio, into=c("target_in", "target_total")) %>% 
    separate(BgRatio, into=c("bg_in", "bg_total")) %>%
    mutate(newpval = ifelse(is.na(pvalue), 1, pvalue), 
           logp = -log10(newpval),
           target_frac = as.numeric(target_in)/as.numeric(target_total), 
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           enrich = log2(target_frac/bg_frac)) %>%
    dplyr::select(version, set, ID, logp, target_frac, bg_frac, enrich) %>% 
    filter(target_frac < 1, bg_frac < 1) %>%
    arrange(ID)
}


orig_bp <- read_tsv(paste0(path_to_data, "/original_GO_galgal4_run1_BP_real.tsv.gz")) %>% transform_real()
orig_mf <- read_tsv(paste0(path_to_data, "/original_GO_galgal4_run1_MF_real.tsv.gz")) %>% transform_real()

#read perms
orig_bp_perm<-readRDS(paste0(path_to_data, "/goperms/original_galgal4_BP.robj"))
orig_mf_perm<-readRDS(paste0(path_to_data, "/goperms/original_galgal4_MF.robj"))

orig_bp_merge <- full_join(orig_bp, orig_bp_perm, by=c("version" = "version", "set" = "set", "ID" = "ID")) 
orig_mf_merge <- full_join(orig_mf, orig_mf_perm, by=c("version" = "version", "set" = "set", "ID" = "ID"))

#compute P-values

orig_bp_merge <-
  orig_bp_merge %>% 
  rowwise %>% 
  mutate(pval_frac = max(1-ecdf_frac(target_frac), 0.0002), 
         pval_logp = max(1-ecdf_logp(logp), 0.0002), 
         pval_enrich = max(1-ecdf_enrich(enrich), 0.0002)) %>% 
  ungroup %>% group_by(version, set) %>% 
  mutate(qval_logp = p.adjust(pval_logp, "BH"),
         qval_frac = p.adjust(pval_frac, "BH"),
         qval_enrich = p.adjust(pval_enrich, "BH"))

orig_mf_merge <-
  orig_mf_merge %>% 
  rowwise %>% 
  mutate(pval_frac = max(1-ecdf_frac(target_frac), 0.0002), 
         pval_logp = max(1-ecdf_logp(logp), 0.0002), 
         pval_enrich = max(1-ecdf_enrich(enrich), 0.0002)) %>% 
  ungroup %>% group_by(version, set) %>% 
  mutate(qval_logp = p.adjust(pval_logp, "BH"),
         qval_frac = p.adjust(pval_frac, "BH"),
         qval_enrich = p.adjust(pval_enrich, "BH"))

#analysis

orig_mf_merge %>% filter(version == "gain") %>% filter(qval_frac < 0.25) %>% 
  mutate(exp_frac = map(ecdf_frac, summary) %>% map_dbl(., 3)) %>%
  select(set, ID, target_frac, exp_frac, bg_frac, enrich, pval_frac, qval_frac) %>% 
  write_tsv(paste0(path_to_data, "/GOPERM_orig_gain_mf_results.tsv"))

orig_bp_merge %>%filter(version == "gain") %>% filter(qval_frac < 0.25) %>% 
  mutate(exp_frac = map(ecdf_frac, summary) %>% map_dbl(., 3)) %>%
  select(set, ID, target_frac, exp_frac, bg_frac, enrich, pval_frac, qval_frac) %>% 
  write_tsv(paste0(path_to_data, "/GOPERM_orig_gain_bp_results.tsv"))
