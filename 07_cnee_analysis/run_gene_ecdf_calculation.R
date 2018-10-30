#revised Oct 2018 for paper revisions
library(tidyverse)

#loop to compute ecdfs for each run (extended, original, reduced, extended cormorant and version galgal4, galgal5)

wdir<-getwd()

for (whichset in c("extended", "original", "reduced", "extended_ratiteVcorm")) {
  for (whichgenome in c("galgal4", "galgal5")) {
    spec_patt<-glob2rx(paste0(whichset, "_gene_", whichgenome, "_run*_perm.tsv"))
    files<-list.files(path=paste0(wdir, "/geneperms"), pattern=spec_patt, full.names = TRUE)
    results<-list()
    for (file in files) {
      results[[file]] <- read_tsv(file) 
    }
    bind_rows(results) %>% group_by(version, set, gene) %>% summarize(ecdf_gene = list(ecdf(rand_TRUE))) %>% saveRDS(file=paste0("geneperms/", whichset, "_", whichgenome, ".robj"))
  }
}