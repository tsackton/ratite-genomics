setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/")
library(tidyverse)

cnee <- read_tsv("cnees.tsv")

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
