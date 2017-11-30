setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/")

obs_bp <- read_tsv("~/Projects/birds/ratite_compgen/data/cnee_perms/obs_bp_results.tsv")
perm_bp <-read_tsv("~/Projects/birds/ratite_compgen/data/cnee_perms/perm_bp_results.tsv")
obs_mf <- read_tsv("~/Projects/birds/ratite_compgen/data/cnee_perms/obs_mf_results.tsv")
perm_mf <-read_tsv("~/Projects/birds/ratite_compgen/data/cnee_perms/perm_mf_results.tsv")

#add stats to obs
obs_bp <- obs_bp %>% mutate(newpval = ifelse(is.na(pvalue), 1, pvalue), logp.perm = -1*log10(newpval)) %>% 
  separate(GeneRatio, into=c("target_in", "target_total")) %>% 
  separate(BgRatio, into=c("bg_in", "bg_total")) %>%
  mutate(target_frac = as.numeric(target_in)/as.numeric(target_total), bg_frac = as.numeric(bg_in)/as.numeric(bg_total))
  
obs_mf <- obs_mf %>% mutate(newpval = ifelse(is.na(pvalue), 1, pvalue), logp.perm = -1*log10(newpval)) %>% 
  separate(GeneRatio, into=c("target_in", "target_total")) %>% 
  separate(BgRatio, into=c("bg_in", "bg_total")) %>%
  mutate(target_frac = as.numeric(target_in)/as.numeric(target_total), bg_frac = as.numeric(bg_in)/as.numeric(bg_total))

#compute p-values for each set in obs bp
get_empirical_pval <- function(term, set, permDF, value, column) {
  val_to_get <- enquo(column)
  null <- permDF %>% filter(set == set, ID == term) %>% pull(!!val_to_get)
  (sum(null >= value)+1) / length(null)
}

#generate
obs_bp_real_logp <- obs_bp %>% group_by(set, ID) %>% mutate(epval = get_empirical_pval(ID, set, perm_bp, logp.perm, logp.perm))
write_tsv(obs_bp_real_logp, path="obs_bp_real_logp.results")
obs_bp_real_targetfrac <- obs_bp %>% group_by(set, ID) %>% mutate(epval = get_empirical_pval(ID, set, perm_bp, target_frac, target_frac))
write_tsv(obs_bp_real_targetfrac, path="obs_bp_real_targetfrac.results")
obs_mf_real_logp <- obs_mf %>% group_by(set, ID) %>% mutate(epval = get_empirical_pval(ID, set, perm_mf, logp.perm, logp.perm))
write_tsv(obs_mf_real_logp, path="obs_mf_real_logp.results")
obs_mf_real_targetfrac <- obs_mf %>% group_by(set, ID) %>% mutate(epval = get_empirical_pval(ID, set, perm_mf, target_frac, target_frac))
write_tsv(obs_mf_real_targetfrac, path="obs_mf_real_targetfrac.results")

#read
obs_bp_real_logp <- read_tsv("obs_bp_real_logp.results")
obs_bp_real_targetfrac <- read_tsv("obs_bp_real_targetfrac.results")
obs_mf_real_logp <- read_tsv("obs_mf_real_logp.results")
obs_mf_real_targetfrac <- read_tsv("obs_mf_real_targetfrac.results")

obs_bp_real_logp <- obs_bp_real_logp %>% ungroup %>% group_by(set) %>% mutate(eqval = p.adjust(epval, "fdr"))
obs_bp_real_targetfrac <- obs_bp_real_targetfrac %>% ungroup %>% group_by(set) %>% mutate(eqval = p.adjust(epval, "fdr"))
obs_mf_real_logp <- obs_mf_real_logp %>% ungroup %>% group_by(set) %>% mutate(eqval = p.adjust(epval, "fdr"))
obs_mf_real_targetfrac <- obs_mf_real_targetfrac %>% ungroup %>% group_by(set) %>% mutate(eqval = p.adjust(epval, "fdr"))

obs_bp_real_logp %>% ungroup %>% filter(eqval < 0.2)
obs_mf_real_logp %>% ungroup %>% filter(eqval < 0.2)
obs_bp_real_targetfrac %>% ungroup %>% filter(eqval < 0.2, qvalue < 0.01) %>% select(set, Description, ID, target_frac, epval, eqval)
obs_mf_real_targetfrac %>% ungroup %>% filter(eqval < 0.2, qvalue < 0.01) %>% select(set, Description, ID, qvalue, target_frac, bg_frac, epval, eqval)

#main go results: DNA binding / nucleic acid binding overrepresented
perm_mf %>% filter(ID=="GO:0003677", set=="set3") %>% summarize(ave = mean(target_frac))
perm_mf %>% filter(ID=="GO:0003676", set=="set3") %>% summarize(ave = mean(target_frac))
perm_bp %>% filter(ID=="GO:0010468", set=="set3") %>% summarize(ave = mean(target_frac))
perm_bp %>% filter(ID=="GO:1903506", set=="set3") %>% summarize(ave = mean(target_frac))


