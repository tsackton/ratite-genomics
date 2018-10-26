library(tidyverse)
library(parallel)

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

#function to do work

compute_gene_results <- function(DF, outname, CORES, PERMS) {
  
  gene_perm_all <- list()
  gene_res_all <- list()
  
  for (ver in c("gain", "gap", "gain_gap", "orig")) {
    
    cnee <- DF %>% filter(version == ver)
    gene_res_all[[ver]] <- lapply(c(quo(rar),quo(crar),quo(crar_dollo)), get_gene_counts, DF=cnee) %>% bind_rows(.id="set")
  
    para_cores <- CORES
    num_perms <- PERMS
  
    gene_counts_perm_list <- list(rar = mclapply(1:num_perms, perm_gene_counts, DF=cnee, indicator=quo(rar), mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"),
                                crar = mclapply(1:num_perms, perm_gene_counts, DF=cnee, indicator=quo(crar), mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"),
                                crar_dollo = mclapply(1:num_perms, perm_gene_counts, DF=cnee, indicator=quo(crar_dollo), mc.preschedule = FALSE, mc.cores = para_cores) %>% bind_rows(.id="perm"))
  
    gene_perm_all[[ver]] <- gene_counts_perm_list %>% bind_rows(.id="set")
  
  }

  bind_rows(gene_perm_all, .id="version") %>% write_tsv(paste0(outname, "_perm.tsv"))
  bind_rows(gene_res_all, .id="version") %>% write_tsv(paste0(outname, "_real.tsv"))
                                                                                                              
}

args <- commandArgs(trailingOnly = TRUE)
#args are 1 annotation file name, 2 permutation index ID for slurm batch processing, 3 number of cores, 4 number of permutations, 

gene_gg<-read_tsv(paste0("../04_wga/03_ce_annotation/cnees.", args[1], ".annotation"), col_names = c("cnee", "gene"))

cnee_orig <- read_tsv("final_original_cnee.tsv.gz") %>% 
  select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
  full_join(gene_gg, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
         crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
  distinct(cnee, version, .keep_all=TRUE) %>%
  select(cnee, version, rar, crar, crar_dollo, gene)

cnee_red <- read_tsv("final_reduced_cnee.tsv.gz") %>% 
  select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
  full_join(gene_gg, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
         crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
  distinct(cnee, version, .keep_all=TRUE) %>%
  select(cnee, version, rar, crar, crar_dollo, gene)

cnee_ext <- read_tsv("final_extended_cnee.tsv.gz") %>% 
  select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
  full_join(gene_gg, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
         crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
  distinct(cnee, version, .keep_all=TRUE) %>%
  select(cnee, version, rar, crar, crar_dollo, gene)

#note in ext2, convergence defined as ratites + cormorants
cnee_ext2 <- read_tsv("final_extended_cnee.tsv.gz") %>% 
  select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo, gc_pp_loss) %>%
  full_join(gene_gg, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
         crar = ifelse(rar & floss_cl_pp >= 1.8 & gc_pp_loss > 0.90, TRUE, FALSE),
         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8 & gc_pp_loss > 0.90, TRUE, FALSE)) %>%
  distinct(cnee, version, .keep_all=TRUE) %>%
  select(cnee, version, rar, crar, crar_dollo, gene)

compute_gene_results(cnee_orig, paste0("geneperms/original_gene_", args[1], "_run", args[2]), args[3], args[4])
compute_gene_results(cnee_ext, paste0("geneperms/extended_gene_", args[1], "_run", args[2]), args[3], args[4])
compute_gene_results(cnee_ext2, paste0("geneperms/extended_ratiteVcorm_gene_", args[1], "_run", args[2]), args[3], args[4])
compute_gene_results(cnee_red, paste0("geneperms/reduced_gene_", args[1], "_run", args[2]), args[3], args[4])


