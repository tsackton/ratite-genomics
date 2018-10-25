### GENE ENRICHMENT HERE ###
#for gene enrichment tests, the idea is to randomly permute each set and get count of CNEEs per gene
#strategy here is to make an indicator variable, compute "real" T/F per gene, and then shuffle indicator variable

get_gene_counts <- function(DF, indicator) {
  #  indcol <- enquo(indicator)
  #  indcol <- indicator
  DF %>% filter(gene != ".") %>% mutate(in_target = !!indicator) %>% 
    count(gene, in_target) %>% filter(!is.na(in_target)) %>% spread(in_target, n, fill=0, drop=FALSE, sep="_")
}

perm_gene_counts <- function(perm, DF, indicator) {
  #  indcol <- enquo(indicator)
  DF %>% filter(gene != ".") %>% mutate(rand = sample(!!indicator)) %>% count(gene, rand) %>% filter(!is.na(rand)) %>% spread(rand, n, fill=0, drop=FALSE, sep="_")
}

gene_gg4<-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/04_wga/03_ce_annotation/cnees.galgal4.annotation", col_names = c("cnee", "gene"))

cnee_orig <- read_tsv("final_original_cnee.tsv.gz") %>% 
  select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
  full_join(gene_gg4, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
         crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
  distinct(cnee, version, .keep_all=TRUE) %>%
  select(cnee, version, rar, crar, crar_dollo, gene)

cnee_red <- read_tsv("final_reduced_cnee.tsv.gz") %>% 
  select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
  full_join(gene_gg4, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
         crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
  distinct(cnee, version, .keep_all=TRUE) %>%
  select(cnee, version, rar, crar, crar_dollo, gene)

cnee_ext <- read_tsv("final_extended_cnee.tsv.gz") %>% 
  select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
  full_join(pos_gg4_uscs, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
         crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
  distinct(cnee, version, .keep_all=TRUE) %>%
  select(cnee, version, rar, crar, crar_dollo, gene)

#note in ext2, convergence defined as ratites + cormorants
cnee_ext2 <- read_tsv("final_extended_cnee.tsv.gz") %>% 
  select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo, gc_pp_loss) %>%
  full_join(gene_gg4, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
         crar = ifelse(rar & floss_cl_pp >= 1.8 & gc_pp_loss > 0.90, TRUE, FALSE),
         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8 & gc_pp_loss > 0.90, TRUE, FALSE)) %>%
  distinct(cnee, version, .keep_all=TRUE) %>%
  select(cnee, version, rar, crar, crar_dollo, gene)


#make real dataset-- this needs editing

gene_counts_obs <- lapply(c(quo(set1),quo(set2),quo(set3),quo(set4)), get_gene_counts, DF=cnee) %>% bind_rows(.id="set")

para_cores <- CORES
num_perms <- PERMS

gene_counts_perm_list <- list(set1 = mclapply(1:num_perms, perm_gene_counts, DF=cnee, indicator=quo(set1), mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"),
                              set2 = mclapply(1:num_perms, perm_gene_counts, DF=cnee, indicator=quo(set2), mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"),
                              set3 = mclapply(1:num_perms, perm_gene_counts, DF=cnee, indicator=quo(set3), mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"),
                              set4 = mclapply(1:num_perms, perm_gene_counts, DF=cnee, indicator=quo(set4), mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"))

gene_counts_perm <- gene_counts_perm_list %>% bind_rows(.id="set")

write_tsv(gene_counts_perm, path="perm_gene_count_results.tsv")
write_tsv(gene_counts_obs, path="obs_gene_count_results.tsv")

#end GENE PERMUTATIONS#