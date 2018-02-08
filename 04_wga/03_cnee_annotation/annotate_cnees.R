library(dplyr)
library(tidyr)
setwd("~/Projects/birds/ratite_compgen/ratite-genomics/annotation/cnee/")
galgal4<-read.table("GCF_000002315.3_Gallus_gallus-4.0_feature_table.txt", stringsAsFactors = FALSE, comment.char = "", header = TRUE, sep="\t", quote="") %>% tbl_df
galgal5<-read.table("GCF_000002315.4_Gallus_gallus-5.0_feature_table.txt", stringsAsFactors = FALSE, comment.char = "", header = TRUE, sep="\t", quote="") %>% tbl_df

galgal4tss <- galgal4 %>% rename(feature = X..feature) %>% 
  filter(feature == "gene", class == "protein_coding") %>% 
  select(genomic_accession, start, end, strand, symbol, GeneID, name) %>%
  mutate(pos = ifelse(strand == "+", start-1, end-1), pos2 = pos+1) %>%
  unite(name2, GeneID, symbol, sep=":") %>%
  select(genomic_accession, pos, pos2, name2)

galgal5tss <- galgal5 %>% rename(feature = X..feature) %>% 
  filter(feature == "gene", class == "protein_coding") %>% 
  select(genomic_accession, start, end, strand, symbol, GeneID, name) %>%
  mutate(pos = ifelse(strand == "+", start-1, end-1), pos2 = pos+1) %>%
  unite(name2, GeneID, symbol, sep=":") %>%
  select(genomic_accession, pos, pos2, name2)

write.table(galgal4tss, "galgal4_tss.bed", sep="\t", col.names = FALSE, row.names = FALSE, quote=F)
write.table(galgal5tss, "galgal5_tss.bed", sep="\t", col.names = FALSE, row.names = FALSE, quote=F)

