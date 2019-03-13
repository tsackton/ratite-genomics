#revised Oct 2018 for paper revisions
library(tidyverse)

#loop to compute ecdfs for each run (extended, original, reduced, extended cormorant and version galgal4, galgal5)

wdir<-"/n/holylfs/LABS/edwards_lab/tsackton/RATITE_PAPER_DATA_FREEZE/DRYAD/07_cnees/processed"

for (whichset in c("extended", "original", "reduced")) {
  for (whichgenome in c("galgal4", "galgal5")) {
    for (whichont in c("BP", "MF")) {
    	spec_patt<-glob2rx(paste0(whichset, "_GO_", whichgenome, "_run*_", whichont, "_perm.tsv"))
    	files<-list.files(path=paste0(wdir, "/goperms"), pattern=spec_patt, full.names = TRUE)
    	results<-list()
    	for (file in files) {
    	  results[[file]] <- read_tsv(file) %>% filter(!is.na(target_frac), target_frac < 1) %>% mutate(enrich = log2(target_frac/bg_frac))
    	}
    	bind_rows(results) %>% 
    	  group_by(version, set, ID) %>% 
    	  summarize(ecdf_frac = list(ecdf(target_frac)),
    	            ecdf_enrich = list(ecdf(enrich)),
    	            ecdf_logp = list(ecdf(logp.perm))) %>% 
    	  saveRDS(file=paste0(wdir, "/goperms/", whichset, "_", whichgenome, "_", whichont, ".robj"))
    	}
  }
}
