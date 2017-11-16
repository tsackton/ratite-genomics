#get qc stats from hog matrix

library(tidyverse)
library(rlang)

hogs<-read_delim("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/all_hog_info.tsv", delim="\t")
full_matrix<-read_delim("~/Projects/birds/ratite_compgen/ratite-genomics/homology/updated_hog_matrix.txt", delim="\t")

sum(full_matrix$aptHaa)
sum(full_matrix$aptRow)
sum(full_matrix$aptOwe)
sum(full_matrix$droNov)
sum(full_matrix$casCas)
sum(full_matrix$notPer)
sum(full_matrix$eudEle)
sum(full_matrix$cryCin)
sum(full_matrix$rheAme)
sum(full_matrix$rhePen)

get_orth_stats <- function(df, tarSp) {
  tarVar <- sym(tarSp)
  targetTotal <- df %>% select(!!tarVar) %>% sum
  scGenes <- df %>% filter(UQ(tarVar) == 1, galGal == 1) %>% select(!!tarVar) %>% sum
  allGenes <- df %>% filter(UQ(tarVar) >= 1, galGal >= 1) %>% select(!!tarVar) %>% sum
  return(data.frame(sp = tarSp, sc=scGenes/targetTotal, all=allGenes / targetTotal))
}

for (name in colnames(full_matrix)[2:43]) {
  print(get_orth_stats(full_matrix, name))
}


