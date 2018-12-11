library(seqinr)
source("~/Projects/genomics/phyloAcc/PhyloAcc/R/drawAlign_function.R")

### Read in tree data
treeData <- prepare_data(tree_path = "~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/neut_original.named.mod", species_name = "~/Projects/genomics/phyloAcc/PhyloAcc/Data/ratite/species_names.txt", common_name = "~/Projects/genomics/phyloAcc/PhyloAcc/Data/ratite/birdname2.txt")

### Generate evolutionary pattern and sequence alignment for one element from PhyloAcc outputs 
#### read in BF scores as well as marginal log likelihood under null, accelerated and full model #### 
score <- read.table("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/03_phyloAcc_res/ORG_1012/gain/Combined_elem_lik.txt", header=T)
## compute BF2
score$BF2 <- score$logBF2
score$log_ratio <- score$logBF1
## order score by BF1
score <- score[order(-score$log_ratio),]

#### read in posteriors of substitution rates and latent conservation states ####
postZ <- read.table("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/03_phyloAcc_res/ORG_1012/gain/Combined_post_Z_M2.txt", header=T, check.names = F)
postZ <- as_tibble(postZ, rownames = "id")

sel <- "mCE967994"
sel <- "mCE1140641"
lk = score %>% filter(ID == sel)
targets = c("strCam","rhePen","rheAme","casCas","droNov","aptRow","aptHaa","aptOwe","anoDid") # target species
Z = postZ %>% filter(id == sel) %>% select(-id) %>% unlist

tit = paste("cnee: ", sel, "logBF1:", round(lk$log_ratio), "logBF2:",round(lk$BF2), "  ") # use BF scores and posterior substitution rates as title
plotZPost(Z, treeData, target_species=targets, tit=tit, offset=5,cex.score = 2) # offset= 6 indicates the posterior of Z start from 7th column

