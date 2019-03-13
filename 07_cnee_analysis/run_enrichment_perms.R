#this code runs the permutations to test for GO, gene, and spatial enrichment in RARs and cRARs
library(tidyverse)
library(clusterProfiler)
library(org.Gg.eg.db)
library(parallel)
library(rlist)

#set working directory
#setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/")

compute_go_results <- function(DF, outname, CORES, PERMS) {

  ##INTERNAL FUNCTIONS##
  calc_enrich <- function(targetset, background, ont) { 
    enrichGO(targetset$ncbi,'org.Gg.eg.db',
             pvalueCutoff=1.5,
             qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, 
             pAdjustMethod="none",
             universe=background,
             keyType="ENTREZID",
             ont=ont) 
  }
  
  get_go_perm <- function(DF, samples, golist, ont) {
    rand <- DF %>% 
      sample_n(samples) %>% 
      filter(gene != ".") %>% 
      dplyr::select(gene) %>% 
      separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
      distinct(ncbi)
    
    
    background <- DF %>% 
      filter(gene != ".") %>% 
      dplyr::select(gene) %>% 
      separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
      distinct(ncbi)
    
    rand_go <- calc_enrich(targetset=rand, background=background$ncbi, ont=ont)
    
    golist %>% 
      left_join(rand_go@result, by=c("ID" = "ID")) %>% 
      separate(GeneRatio, into=c("target_in", "target_total")) %>% 
      separate(BgRatio, into=c("bg_in", "bg_total")) %>%
      mutate(newpval = ifelse(is.na(pvalue), 1, pvalue), 
             logp.perm = -log10(newpval),
             target_frac = as.numeric(target_in)/as.numeric(target_total), 
             bg_frac = as.numeric(bg_in)/as.numeric(bg_total)) %>%
      dplyr::select(ID, logp.perm, target_frac, bg_frac) %>% 
      arrange(ID)
  }
  
  #to do one permutation of the full set
  get_one_perm_set <- function(perm, input, DF, golist, ont) {
    try(lapply(input, get_go_perm, DF=DF, golist=golist, ont=ont) %>%
          dplyr::bind_rows(.id="set"), TRUE)
  }
  
  bp_perm_all <- list()
  bp_res_all <- list()

  mf_perm_all <- list()
  mf_res_all <- list()
  
  for (ver in c("gain", "gap", "gain_gap", "orig")) {
  
    cnee <- DF %>% filter(version == ver)
    
    background <- cnee %>% 
      filter(gene != ".") %>% 
      dplyr::select(gene) %>% 
      separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
      distinct(ncbi)
    
    set1 <- cnee %>% filter(gene != ".", rar) %>% 
      dplyr::select(gene) %>% 
      separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
      distinct(ncbi)
    
    set2 <- cnee %>% filter(gene != ".", crar) %>% 
      dplyr::select(gene) %>% 
      separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
      distinct(ncbi)
    
    set3 <- cnee %>% filter(gene != ".", crar_dollo) %>% 
      dplyr::select(gene) %>% 
      separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
      distinct(ncbi)
    
    inputs <- list("rar" = set1, "crar" = set2, "crar_dollo" = set3)
    
      
    bp_res_all[[ver]] <- lapply(inputs, calc_enrich, background=background$ncbi, ont="BP") %>% 
      lapply(slot, name="result") %>% 
      dplyr::bind_rows(.id = "set")
    
    mf_res_all[[ver]] <-  lapply(inputs, calc_enrich, background=background$ncbi, ont="MF") %>% 
      lapply(slot, name="result") %>% 
      dplyr::bind_rows(.id = "set")
    
    merged_mf_terms <- mf_res_all[[ver]] %>% dplyr::distinct(ID)
    merged_bp_terms <- bp_res_all[[ver]] %>% dplyr::distinct(ID)
  
    input_counts<-list("rar" = cnee %>% filter(gene != ".", rar) %>% count %>% pull(n),
                       "crar" = cnee %>% filter(gene != ".", crar) %>% count %>% pull(n),
                       "crar_dollo" = cnee %>% filter(gene != ".", crar_dollo) %>% count %>% pull(n))
    
    bp_perm_all[[ver]] <- mclapply(1:PERMS, get_one_perm_set, input=input_counts, DF=cnee, golist=merged_bp_terms, ont="BP", mc.cores=CORES, mc.preschedule = FALSE) %>%
      list.filter(class(.) == "data.frame") %>% 
      dplyr::bind_rows(.id="perm")
    
    mf_perm_all[[ver]] <- mclapply(1:PERMS, get_one_perm_set, input=input_counts, DF=cnee, golist=merged_mf_terms, ont="MF", mc.cores=CORES, mc.preschedule = FALSE) %>%
      list.filter(class(.) == "data.frame") %>% 
      dplyr::bind_rows(.id="perm")
    
  }

  bind_rows(bp_perm_all, .id="version") %>% write_tsv(paste0(outname, "_BP_perm.tsv"))
  bind_rows(bp_res_all, .id="version") %>% write_tsv(paste0(outname, "_BP_real.tsv"))              
  bind_rows(mf_perm_all, .id="version") %>% write_tsv(paste0(outname, "_MF_perm.tsv"))
  bind_rows(mf_res_all, .id="version") %>% write_tsv(paste0(outname, "_MF_real.tsv"))       
}

### REAL WORK ###

args <- commandArgs(trailingOnly = TRUE)
#args are 1 annotation file name, 2 permutation index ID for slurm batch processing, 3 number of cores, 4 number of permutations, 5 is data path

path_to_data <- args[5]

gene_gg<-read_tsv(paste0("../04_wga/03_ce_annotation/cnees.", args[1], ".annotation"), col_names = c("cnee", "gene"))

cnee_orig <- read_tsv(paste0(path_to_data, "/final_original_cnee.tsv.gz")) %>% 
  dplyr::select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
  full_join(gene_gg, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
         crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
  distinct(cnee, version, .keep_all=TRUE) %>%
  dplyr::select(cnee, version, rar, crar, crar_dollo, gene)

cnee_red <- read_tsv(paste0(path_to_data, "/final_reduced_cnee.tsv.gz")) %>% 
  dplyr::select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
  full_join(gene_gg, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
         crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
  distinct(cnee, version, .keep_all=TRUE) %>%
  dplyr::select(cnee, version, rar, crar, crar_dollo, gene)

cnee_ext <- read_tsv(paste0(path_to_data, "/final_extended_cnee.tsv.gz")) %>% 
  dplyr::select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
  full_join(gene_gg, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
         crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
  distinct(cnee, version, .keep_all=TRUE) %>%
  dplyr::select(cnee, version, rar, crar, crar_dollo, gene)

#note in ext2, convergence defined as ratites + cormorants
#cnee_ext2 <- read_tsv("final_extended_cnee.tsv.gz") %>% 
#  dplyr::select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo, gc_pp_loss) %>%
#  full_join(gene_gg, by=c("cnee" = "cnee")) %>%
#  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
#         crar = ifelse(rar & floss_cl_pp >= 1.8 & gc_pp_loss > 0.90, TRUE, FALSE),
#         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8 & gc_pp_loss > 0.90, TRUE, FALSE)) %>%
#  distinct(cnee, version, .keep_all=TRUE) %>%
#  dplyr::select(cnee, version, rar, crar, crar_dollo, gene)

compute_go_results(cnee_orig, paste0(path_to_data, "/goperms/original_GO_", args[1], "_run", args[2]), args[3], args[4])
compute_go_results(cnee_ext, paste0(path_to_data, "/goperms/extended_GO_", args[1], "_run", args[2]), args[3], args[4])
#compute_go_results(cnee_ext2, paste0(path_to_data, "/goperms/extended_ratiteVcorm_GO_", args[1], "_run", args[2]), args[3], args[4])
compute_go_results(cnee_red, paste0(path_to_data, "/goperms/reduced_GO_", args[1], "_run", args[2]), args[3], args[4])
