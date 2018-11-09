#make transfac datasets
library(tidyverse)
setwd("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/transfac_input/")

scores <- read_tsv("cnee_gapmissing_extended_ncbi_galgal4.bed", col_names=c("chr", "start", "end", "id", "score")) %>% select(id, score)

part_bed <- read_tsv("allspecies_cnee_concat_partitions.bed", col_names = c("chr", "start", "end", "id"))

part_bed %>% full_join(scores, by = c("id" = "id")) %>% 
  mutate(chr = "SPECIES") %>%
  filter(score != "none") %>% select(chr, start, end) %>%
  write_tsv(path="cnee_part_acc.bed", col_names = FALSE)
