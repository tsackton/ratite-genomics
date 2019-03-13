#revised Nov 2018 for paper revisions
library(tidyverse)

transform_real <- function(DF) {
  DF %>% 
    separate(GeneRatio, into=c("target_in", "target_total")) %>% 
    separate(BgRatio, into=c("bg_in", "bg_total")) %>%
    mutate(newpval = ifelse(is.na(pvalue), 1, pvalue), 
           logp = -log10(newpval),
           target_frac = as.numeric(target_in)/as.numeric(target_total), 
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total)) %>%
    dplyr::select(version, set, ID, logp, target_frac, bg_frac) %>% 
    arrange(ID)
}


orig_bp <- read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/goperms/original_GO_galgal4_run1_BP_real.tsv") %>% transform_real()

orig_mf <- read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/goperms/original_GO_galgal4_run1_MF_real.tsv") %>% transform_real()

#read perms
orig_bp_perm<-readRDS("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/goperms/original_galgal4_BP.robj")
orig_mf_perm<-readRDS("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/goperms/original_galgal4_MF.robj")

#compute P-values

orig_bp_merge <- full_join(orig_bp, orig_bp_perm, by=c("version" = "version", "set" = "set", "ID" = "ID")) %>% 
  rowwise %>% 
  mutate(pval_frac = max(1-ecdf_frac(target_frac), 1e-04), 
         pval_logp = max(1-ecdf_logp(logp), 1e-04), 
         pval_enrich = max(1-ecdf_enrich(log2(target_frac/bg_frac))), 1e-04) %>% 
  ungroup %>% group_by(version, set) %>% 
  mutate(qval_logp = p.adjust(pval_logp, "BH"),
         qval_frac = p.adjust(pval_frac, "BH"),
         qval_enrich = p.adjust(pval_enrich, "BH"))


orig_mf_merge <- full_join(orig_mf, orig_mf_perm, by=c("version" = "version", "set" = "set", "ID" = "ID")) %>% 
  rowwise %>% 
  mutate(pval_frac = max(1-ecdf_frac(target_frac), 1e-04), 
         pval_logp = max(1-ecdf_logp(logp), 1e-04), 
         pval_enrich = max(1-ecdf_enrich(log2(target_frac/bg_frac)), 1e-04)) %>% 
  ungroup %>% group_by(version, set) %>% 
  mutate(qval_logp = p.adjust(pval_logp, "BH"),
         qval_frac = p.adjust(pval_frac, "BH"),
         qval_enrich = p.adjust(pval_enrich, "BH"))

orig_mf_merge %>% filter(version == "gain") %>% filter(qval_frac < 0.20) %>% 
  select(set, ID, target_frac, bg_frac, pval_frac, qval_frac) %>% 
  write_tsv("~/Projects/birds/ratite_compgen/manuscript/ScienceSubmissionRev2/orig_gain_mf_results_2.tsv")
  
#ugly hack
orig_mf_merge %>% filter(version == "gain") %>% filter(qval_frac < 0.20) %>% pull(ecdf_frac) %>% map(., summary) %>% map_dbl(., 3)
  
  
orig_bp_merge %>% filter(version == "gain") %>% filter(qval_frac < 0.20)  %>% 
  select(set, ID, target_frac, bg_frac, pval_frac, qval_frac)  %>% 
  write_tsv("~/Projects/birds/ratite_compgen/manuscript/ScienceSubmissionRev2/orig_gain_bp_results_2.tsv")

#ugly hack
orig_bp_merge %>% filter(version == "gain") %>% filter(qval_frac < 0.20) %>% pull(ecdf_frac) %>% map(., summary) %>% map_dbl(., 3) %>% as.data.frame() %>%
  write_tsv("~/Projects/birds/ratite_compgen/manuscript/ScienceSubmissionRev2/orig_gain_bp_results_2-etf.tsv")



summary((orig_mf_merge %>% filter(ID == "GO:0003676", version == "gain", set == "crar") %>% pull(ecdf_frac))[[1]])
summary((orig_mf_merge %>% filter(ID == "GO:0003677", version == "gain", set == "crar") %>% pull(ecdf_frac))[[1]])
summary((orig_mf_merge %>% filter(ID == "GO:0043565", version == "gain", set == "crar") %>% pull(ecdf_frac))[[1]])
summary((orig_bp_merge %>% filter(ID == "GO:0006355", version == "gain", set == "crar") %>% pull(ecdf_frac))[[1]])
summary((orig_bp_merge %>% filter(ID == "GO:0060173", version == "gain", set == "crar") %>% pull(ecdf_frac))[[1]])

