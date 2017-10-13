#chick
library(tidyverse)
colNames = c("pattern", "target", "start", "stop", "strand", "score", "pval", "qval", "matched_seq")
chick<-read_delim(file="~/Projects/birds/ratite_compgen/ratite-genomics/analysis/motifs/chick_out/fimo.txt", delim="\t", col_names = colNames, comment="#")
tin<-read_delim(file="~/Projects/birds/ratite_compgen/ratite-genomics/analysis/motifs/tinamou_out/fimo.txt", delim="\t", col_names = colNames, comment="#")
rhea<-read_delim(file="~/Projects/birds/ratite_compgen/ratite-genomics/analysis/motifs/rhea_out/fimo.txt", delim="\t", col_names = colNames, comment="#")
moa<-read_delim(file="~/Projects/birds/ratite_compgen/ratite-genomics/analysis/motifs/moa_out/fimo.txt", delim="\t", col_names = colNames, comment="#")



chick %>% filter(pattern %in% tin$pattern, !(pattern %in% rhea$pattern), !(pattern %in% moa$pattern), qval < 0.2) %>% distinct(pattern) %>% arrange(pattern)
tin %>% filter(pattern %in% chick$pattern, !(pattern %in% rhea$pattern),  !(pattern %in% moa$pattern), qval < 0.2) %>% distinct(pattern) %>% arrange(pattern)

chick %>% filter(pattern %in% rhea$pattern, !(pattern %in% tin$pattern), qval < 0.05) %>% distinct(pattern) %>% arrange(pattern)
rhea %>% filter(pattern %in% chick$pattern, !(pattern %in% tin$pattern), qval < 0.05) %>% distinct(pattern) %>% arrange(pattern)

rhea %>% filter(!(pattern %in% chick$pattern), !(pattern %in% tin$pattern), pattern %in% moa$pattern, qval < 0.05) %>% distinct(pattern) %>% arrange(pattern)
