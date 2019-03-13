#code updated for revisions Oct 2018

setwd("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/") #sorry!
library(tidyverse)
library(GenomicRanges)

pos_gg4_uscs<-read_tsv("../04_wga/03_ce_annotation/cnees.galGal4UCSC.bed", col_names = c("chr", "start", "end", "cnee"))

chr_to_remove <- pos_gg4_uscs %>% count(chr) %>% filter(n < 150) %>% pull(chr)

cnee_orig <- read_tsv("final_original_cnee.tsv.gz") %>% 
  select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
  full_join(pos_gg4_uscs, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
         crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
  distinct(cnee, version, .keep_all=TRUE) %>%
  filter(!(chr %in% chr_to_remove)) %>%
  select(chr, start, end, cnee, version, rar, crar, crar_dollo)

cnee_red <- read_tsv("final_reduced_cnee.tsv.gz") %>% 
  select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
  full_join(pos_gg4_uscs, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
         crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
  distinct(cnee, version, .keep_all=TRUE) %>%
  filter(!(chr %in% chr_to_remove))  %>%
  select(chr, start, end, cnee, version, rar, crar, crar_dollo)

cnee_ext <- read_tsv("final_extended_cnee.tsv.gz") %>% 
  select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, floss_cl_pp, floss_cl_pp_dollo) %>%
  full_join(pos_gg4_uscs, by=c("cnee" = "cnee")) %>%
  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
         crar = ifelse(rar & floss_cl_pp >= 1.8, TRUE, FALSE),
         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8, TRUE, FALSE)) %>%
  distinct(cnee, version, .keep_all=TRUE) %>%
  filter(!(chr %in% chr_to_remove)) %>%
  select(chr, start, end, cnee, version, rar, crar, crar_dollo)

#note in ext2, convergence defined as ratites + cormorants
#cnee_ext2 <- read_tsv("final_extended_cnee.tsv.gz") %>% 
#  select(version, cnee, logBF1, logBF2, it_pp_loss, ti_pp_loss, neo_tip_loss, #floss_cl_pp, floss_cl_pp_dollo, gc_pp_loss) %>%
#  full_join(pos_gg4_uscs, by=c("cnee" = "cnee")) %>%
#  mutate(rar = ifelse(logBF1 >= 10 & logBF2 >= 1 & (it_pp_loss + ti_pp_loss) < 1 & neo_tip_loss < 1, TRUE, FALSE),
#         crar = ifelse(rar & floss_cl_pp >= 1.8 & gc_pp_loss > 0.90, TRUE, FALSE),
#         crar_dollo = ifelse(rar & floss_cl_pp_dollo >= 1.8 & gc_pp_loss > 0.90, #TRUE, FALSE)) %>%
#  distinct(cnee, version, .keep_all=TRUE) %>%
#  filter(!(chr %in% chr_to_remove)) %>%
#  select(chr, start, end, cnee, version, rar, crar, crar_dollo)


#we'll use pbinom to compute a simple test of whether, given random distribution of RARs/cRARS, we'd expect X or more in a window
#q = number RARs in window
#size = number CNEEs in window
#prob = total RARs / total CNEEs (in set)
#lower.tail = F will give P-value of q or greater observed given prob

#need to define windows, then get q and size for each window, and compute p-value

#read in seqlengths
seqlen_df<-read.table("galGal4-seqlengths.txt", header=F, stringsAsFactors = F)
seqlen<-seqlen_df$V2
names(seqlen)<-seqlen_df$V1

#function to do analysis

compute_spatial_results <- function(DF, outname) {
  
  spa_res <- list()
  
  for (ver in c("gain", "gain_gap")) {
    
    cnee <- DF %>% filter(version == ver)
    
    #convert to window
    cnee_ranges<-makeGRangesFromDataFrame(cnee, keep.extra.columns = TRUE, ignore.strand=TRUE, seqinfo = seqlen)
    windows_500kb<-tileGenome(seqinfo(cnee_ranges), tilewidth=500000, cut.last.tile.in.chrom = TRUE)
    windows_1000kb<-tileGenome(seqinfo(cnee_ranges), tilewidth=1000000, cut.last.tile.in.chrom = TRUE)
    windows_100kb<-tileGenome(seqinfo(cnee_ranges), tilewidth=100000, cut.last.tile.in.chrom = TRUE)
    
    #1Mb sliding window
    #hacky
    windows_chr<-tileGenome(seqinfo(cnee_ranges), tilewidth=max(seqlengths(cnee_ranges)), cut.last.tile.in.chrom=TRUE)
    windows_1Mb_slide <- slidingWindows(windows_chr, 1e6, 100000) %>% unlist
    
    rar_prob<-sum(cnee_ranges$rar, na.rm=T)/length(cnee_ranges$rar)
    crar_prob<-sum(cnee_ranges$crar, na.rm=T)/length(cnee_ranges$crar)
    crar_dollo_prob<-sum(cnee_ranges$crar_dollo, na.rm=T)/length(cnee_ranges$crar_dollo)

    get_stats_for_window <- function(query_ranges, window_to_test) {
      window_df <- as.data.frame(window_to_test)
      stopifnot(length(window_df$seqnames)==1)
      subsetByOverlaps(query_ranges, window_to_test) %>% 
        as.data.frame %>% 
        summarize(rar_n = sum(rar, na.rm=T), 
                  crar_n = sum(crar, na.rm=T), 
                  crar_dollo_n = sum(crar_dollo, na.rm=T),
                  total_n=n()) %>%
        mutate(rar_pval = pbinom(q=rar_n, size=total_n, prob=rar_prob, lower.tail=FALSE) + dbinom(x=rar_n, size=total_n, prob=rar_prob),
               crar_pval = pbinom(q=crar_n, size=total_n, prob=crar_prob, lower.tail=FALSE) + dbinom(x=crar_n, size=total_n, prob=crar_prob),
               crar_dollo_pval = pbinom(q=crar_dollo_n, size=total_n, prob=crar_dollo_prob, lower.tail=FALSE) + dbinom(x=crar_dollo_n, size=total_n, prob=crar_dollo_prob),
               chr= window_df$seqnames, start=window_df$start, end=window_df$end) %>% 
        select(chr, start, end, rar_n, crar_n, crar_dollo_n, total_n, rar_pval, crar_pval, crar_dollo_pval)
    }
    
    se_500kb<-lapply(1:length(windows_500kb), function(x) get_stats_for_window(cnee_ranges, windows_500kb[x])) %>% dplyr::bind_rows(.id="window") %>% 
      filter(total_n > 0)
    se_500kb$rar_qval <- p.adjust(se_500kb$rar_pval, "fdr")
    se_500kb$crar_qval <- p.adjust(se_500kb$crar_pval, "fdr")
    se_500kb$crar_dollo_qval <- p.adjust(se_500kb$crar_dollo_pval, "fdr")
    se_500kb$window_size = "500kb"
    
    se_100kb<-lapply(1:length(windows_100kb), function(x) get_stats_for_window(cnee_ranges, windows_100kb[x])) %>% dplyr::bind_rows(.id="window") %>% 
      filter(total_n > 0)
    se_100kb$rar_qval <- p.adjust(se_100kb$rar_pval, "fdr")
    se_100kb$crar_qval <- p.adjust(se_100kb$crar_pval, "fdr")
    se_100kb$crar_dollo_qval <- p.adjust(se_100kb$crar_dollo_pval, "fdr")
    se_100kb$window_size = "100kb"
    
    
    se_1000kb<-lapply(1:length(windows_1000kb), function(x) get_stats_for_window(cnee_ranges, windows_1000kb[x])) %>% dplyr::bind_rows(.id="window") %>% 
      filter(total_n > 0)
    se_1000kb$rar_qval <- p.adjust(se_1000kb$rar_pval, "fdr")
    se_1000kb$crar_qval <- p.adjust(se_1000kb$crar_pval, "fdr")
    se_1000kb$crar_dollo_qval <- p.adjust(se_1000kb$crar_dollo_pval, "fdr")
    se_1000kb$window_size = "1000kb"
    
    #sliding windows
    se_slide<-lapply(1:length(windows_1Mb_slide), function(x) get_stats_for_window(cnee_ranges, windows_1Mb_slide[x])) %>% dplyr::bind_rows(.id="window") %>% 
      filter(total_n > 0)
    se_slide$rar_qval <- p.adjust(se_slide$rar_pval, "fdr")
    se_slide$crar_qval <- p.adjust(se_slide$crar_pval, "fdr")
    se_slide$crar_dollo_qval <- p.adjust(se_slide$crar_dollo_pval, "fdr")
    se_slide$window_size = "1000kb_100kb_slide"
    
    spa_res[[ver]] <- rbind(se_100kb, se_500kb, se_1000kb, se_slide) %>% as.tibble
  }

  bind_rows(spa_res, .id="version") %>% write_tsv(outname)
}

compute_spatial_results(cnee_orig, "original_spatial_results.tsv")
compute_spatial_results(cnee_ext, "extended_spatial_results.tsv")
compute_spatial_results(cnee_red, "reduced_spatial_results.tsv")

