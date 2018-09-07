#this needs some work...?

setwd("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/")
library(tidyverse)

cnee <- read_tsv("cnees.tsv.gz")

#filter

cnee <- select(cnee, chr, start, end, cnee, ratite_accel.1:ratite_conv.2, ratite_loss.mat:ratite_loss_cons_min.mat, ratite_loss.prob:ratite_loss_cons_min.prob) %>%
  distinct(cnee, .keep_all=TRUE)

#remove non-major chromosomes as pointless to test clustering (those with < 500 CNEEs)
cnee <- filter(cnee, !(chr %in% c("chr25","chr1_AADN03009345_random", "chr1_AADN03009420_random", "chr1_AADN03009258_random", "chr1_AADN03009284_random", "chr1_AADN03009295_random", "chr1_AADN03009295_random", "chr1_AADN03009295_random", "chr16", "chrLGE22C19W28_E50C23", "chrW", "chrLGE64")))
               
#we'll use pbinom to compute a simple test of whether, given random distribution of RARs/cRARS, we'd expect X or more in a window
#q = number RARs in window
#size = number CNEEs in window
#prob = total RARs / total CNEEs (in set)
#lower.tail = F will give P-value of q or greater observed given prob

#need to define windows, then get q and size for each window, and compute p-value
library("GenomicRanges")

#read in seqlengths
seqlen_df<-read.table("galGal4-seqlengths.txt", header=F, stringsAsFactors = F)
seqlen<-seqlen_df$V2
names(seqlen)<-seqlen_df$V1

#convert CNEEs to GRanges object
cnee_ranges<-makeGRangesFromDataFrame(cnee, keep.extra.columns = TRUE, ignore.strand=TRUE, seqinfo = seqlen)
windows_500kb<-tileGenome(seqinfo(cnee_ranges), tilewidth=500000, cut.last.tile.in.chrom = TRUE)
windows_1000kb<-tileGenome(seqinfo(cnee_ranges), tilewidth=1000000, cut.last.tile.in.chrom = TRUE)
windows_100kb<-tileGenome(seqinfo(cnee_ranges), tilewidth=100000, cut.last.tile.in.chrom = TRUE)

#1Mb sliding window
windows_1Mb_slide <- trim(resize(unlist(windows_100kb), width=1e6))

rar.1.prob<-sum(cnee_ranges$ratite_accel.1 & cnee_ranges$ratite_spec.1, na.rm=T)/length(cnee_ranges$ratite_spec.1)
crar.1.prob<-sum(cnee_ranges$ratite_accel.1 & cnee_ranges$ratite_spec.1 & cnee_ranges$ratite_conv.1, na.rm=T)/length(cnee_ranges$ratite_spec.1)
rar.2.prob<-sum(cnee_ranges$ratite_accel.2 & cnee_ranges$ratite_spec.2, na.rm=T)/length(cnee_ranges$ratite_spec.1)
crar.2.prob<-sum(cnee_ranges$ratite_accel.1 & cnee_ranges$ratite_spec.1 & cnee_ranges$ratite_conv.1 & cnee_ranges$ratite_loss_cons_min.mat >= 2, na.rm=T)/length(cnee_ranges$ratite_spec.1)

get_stats_for_window <- function(query_ranges, window_to_test) {
  window_df <- as.data.frame(window_to_test)
  stopifnot(length(window_df$seqnames)==1)
  subsetByOverlaps(query_ranges, window_to_test) %>% 
    as.data.frame %>% 
    summarize(rar.1 = sum(ratite_accel.1 & ratite_spec.1, na.rm=T), rar.2=sum(ratite_accel.2 & ratite_spec.2, na.rm=T), crar.1 = sum(ratite_accel.1 & ratite_spec.1 & ratite_conv.1, na.rm=T), crar.2 = sum(ratite_accel.1 & ratite_spec.1 & ratite_conv.1 & ratite_loss_cons_min.mat >= 2, na.rm=T), total=n()) %>%
    mutate(rar.1.pval = pbinom(q=rar.1, size=total, prob=rar.1.prob, lower.tail=FALSE) + dbinom(x=rar.1, size=total, prob=rar.1.prob)) %>%
    mutate(crar.1.pval = pbinom(q=crar.1, size=total, prob=crar.1.prob, lower.tail=FALSE) + dbinom(x=crar.1, size=total, prob=crar.1.prob)) %>%
    mutate(rar.2.pval = pbinom(q=rar.2, size=total, prob=rar.2.prob, lower.tail=FALSE) + dbinom(x=rar.2, size=total, prob=rar.2.prob)) %>%
    mutate(crar.2.pval = pbinom(q=crar.2, size=total, prob=crar.2.prob, lower.tail=FALSE) + dbinom(x=crar.2, size=total, prob=crar.2.prob)) %>%
    mutate(chr = window_df$seqnames, start=window_df$start, end=window_df$end) %>% 
    select(chr, start, end, rar.1, rar.2, crar.1, crar.2, total, rar.1.pval, crar.1.pval, rar.2.pval, crar.2.pval)
}

se_500kb<-lapply(1:length(windows_500kb), function(x) get_stats_for_window(cnee_ranges, windows_500kb[x])) %>% dplyr::bind_rows(.id="window")
se_500kb <- se_500kb %>% filter(total > 0)
se_500kb$rar.1.qval <- p.adjust(se_500kb$rar.1.pval, "fdr")
se_500kb$crar.1.qval <- p.adjust(se_500kb$crar.1.pval, "fdr")
se_500kb$rar.2.qval <- p.adjust(se_500kb$rar.2.pval, "fdr")
se_500kb$crar.2.qval <- p.adjust(se_500kb$crar.2.pval, "fdr")
se_500kb$window_size = "500kb"

se_100kb<-lapply(1:length(windows_100kb), function(x) get_stats_for_window(cnee_ranges, windows_100kb[x])) %>% dplyr::bind_rows(.id="window")
se_100kb <- se_100kb %>% filter(total > 0)
se_100kb$rar.1.qval <- p.adjust(se_100kb$rar.1.pval, "fdr")
se_100kb$crar.1.qval <- p.adjust(se_100kb$crar.1.pval, "fdr")
se_100kb$rar.2.qval <- p.adjust(se_100kb$rar.2.pval, "fdr")
se_100kb$crar.2.qval <- p.adjust(se_100kb$crar.2.pval, "fdr")
se_100kb$window_size = "100kb"


se_1000kb<-lapply(1:length(windows_1000kb), function(x) get_stats_for_window(cnee_ranges, windows_1000kb[x])) %>% dplyr::bind_rows(.id="window")
se_1000kb <- se_1000kb %>% filter(total > 0)
se_1000kb$rar.1.qval <- p.adjust(se_1000kb$rar.1.pval, "fdr")
se_1000kb$crar.1.qval <- p.adjust(se_1000kb$crar.1.pval, "fdr")
se_1000kb$rar.2.qval <- p.adjust(se_1000kb$rar.2.pval, "fdr")
se_1000kb$crar.2.qval <- p.adjust(se_1000kb$crar.2.pval, "fdr")
se_1000kb$window_size = "1000kb"

#sliding windows
se_slide<-lapply(1:length(windows_1Mb_slide), function(x) get_stats_for_window(cnee_ranges, windows_1Mb_slide[x])) %>% dplyr::bind_rows(.id="window")
se_slide <- se_slide %>% filter(total > 0)
se_slide$rar.1.qval <- p.adjust(se_slide$rar.1.pval, "fdr")
se_slide$crar.1.qval <- p.adjust(se_slide$crar.1.pval, "fdr")
se_slide$rar.2.qval <- p.adjust(se_slide$rar.2.pval, "fdr")
se_slide$crar.2.qval <- p.adjust(se_slide$crar.2.pval, "fdr")
se_slide$window_size = "1000kb_100kb_slide"

se_results <- rbind(se_100kb, se_500kb, se_1000kb, se_slide) %>% as.tibble
write_tsv(se_results, path="spatial_enrichment.results")

