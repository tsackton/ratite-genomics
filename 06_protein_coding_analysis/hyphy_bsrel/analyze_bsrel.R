setwd("~/Projects/birds/ratite_compgen/ratite-genomics/06_protein_coding_analysis/hyphy_bsrel/")
library(tidyverse)
library(clusterProfiler)
library(biomaRt)

##ANNOTATION -- code copied from https://github.com/ajshultz/avian-immunity/blob/master/ComparativeGenomics/PAML_Analyze/02_GeneID_Annotation.R ##

#Read in HOG IDs
hog_ids <- read.table("../../03_homology/new_hog_list.txt")
colnames(hog_ids) <- c("HOG2_HogID","NCBI_ID","entrezgene","sp")

#Data cleanup, extract chicken gene IDs
hog_ids_gg <- hog_ids %>% tbl_df %>%
  separate(HOG2_HogID,sep = "_",into=c("drop","hog")) %>%
  mutate(entrezgene = as.character(entrezgene), sp = as.character(sp)) %>%
  filter(sp == "galGal") %>%
  dplyr::select(hog,entrezgene)

#Data cleanup, extract zebra finch gene IDs
hog_ids_zf <- hog_ids %>% tbl_df %>%
  separate(HOG2_HogID,sep = "_",into=c("drop","hog")) %>%
  mutate(entrezgene_zf = as.character(entrezgene), sp = as.character(sp)) %>%
  filter(sp == "taeGut") %>%
  dplyr::select(hog,entrezgene_zf)

####Add human NCBI gene IDs from biomaRt
ensembl = useMart("ensembl")

#Use correct datasets - sometimes throws an error for some reason saying the given dataset is not valid. However, it works if you try again once or twice
ensembl_gg = useDataset("ggallus_gene_ensembl",mart=ensembl)
ensembl_hs = useDataset("hsapiens_gene_ensembl",mart=ensembl)
ensembl_zf = useDataset("tguttata_gene_ensembl",mart=ensembl)

#Get human ensembl IDs for homologs
human_ids_gg <- getBM(attributes=c("ensembl_gene_id","external_gene_name","hsapiens_homolog_ensembl_gene","hsapiens_homolog_orthology_confidence"),
                      filters=c("with_hsapiens_homolog"),
                      values=TRUE,
                      mart=ensembl_gg) %>% tbl_df

human_ids_zf <- getBM(attributes=c("ensembl_gene_id","external_gene_name","hsapiens_homolog_ensembl_gene","hsapiens_homolog_orthology_confidence"),
                      filters=c("with_hsapiens_homolog"),
                      values=TRUE,
                      mart=ensembl_zf) %>% tbl_df %>%
  dplyr::rename(ensembl_gene_id_zf = ensembl_gene_id)


#Get ensembl ID to NCBI entrezgene mapping for three different species
gg_ens_to_entrezgene <- getBM(attributes=c("ensembl_gene_id","entrezgene"),
                              mart=ensembl_gg) %>%
  tbl_df %>%
  filter(!is.na(entrezgene)) %>%
  mutate(entrezgene = as.character(entrezgene))

zf_ens_to_entrezgene <- getBM(attributes=c("ensembl_gene_id","entrezgene"),
                              mart=ensembl_zf) %>%
  tbl_df %>%
  filter(!is.na(entrezgene)) %>%
  dplyr::rename(entrezgene_zf=entrezgene,ensembl_gene_id_zf=ensembl_gene_id) %>%
  mutate(entrezgene_zf = as.character(entrezgene_zf))

hs_ens_to_entrezgene <- getBM(attributes=c("ensembl_gene_id","entrezgene"),
                              mart=ensembl_hs) %>%
  tbl_df %>%
  filter(!is.na(entrezgene)) %>%
  dplyr::rename(entrezgene_hs=entrezgene,hsapiens_homolog_ensembl_gene=ensembl_gene_id) %>%
  mutate(entrezgene_hs = as.character(entrezgene_hs))

#Add human IDs to chicken ensembl to entrezgene table

gg_trans_table <- gg_ens_to_entrezgene %>%
  left_join(human_ids_gg,by="ensembl_gene_id") %>%
  left_join(hs_ens_to_entrezgene,by="hsapiens_homolog_ensembl_gene")

zf_trans_table <- zf_ens_to_entrezgene %>%
  left_join(human_ids_zf,by="ensembl_gene_id_zf") %>%
  left_join(hs_ens_to_entrezgene,by="hsapiens_homolog_ensembl_gene")


hog_to_gene <- gg_trans_table %>% full_join(hog_ids_gg) %>% filter(!is.na(hog))

full_list <- read_tsv("../paml_ancrec/paml_M0_parsed.txt.gz", col_types = "cccc?????") %>% 
  filter(!is.na(hog), model == "ancrec") %>%
  mutate(tree = paste0("tree", treenum)) %>%
  dplyr::select(hog, tree, species_tree)

##LOAD DATA## 
bsrel <- read.table("bsrel_res_parsed_ratites_2018-10-15.txt", fill=TRUE, stringsAsFactors = FALSE)
names(bsrel) <- c("class", "tree", "hog", "tsel.s", "nsel.s", "tsel.n", "nsel.n", "tnon", "nnon", "strict_branches", "nom_branches")
bsrel <- bsrel %>% mutate(hog = as.character(hog)) %>% full_join(full_list, by=c("hog" = "hog", "tree" = "tree")) %>% 
  as.tibble() %>% left_join(hog_to_gene, by=c("hog" = "hog"))
  

##ANALYSIS##
ratite_spec_genes <- bsrel %>% filter(tsel.n >= 1, nsel.n == 0, species_tree == "True") %>% filter(!is.na(ensembl_gene_id))
ratite_spec_genes_gt <- bsrel %>% filter(tsel.n >= 1, nsel.n == 0, species_tree == "False") %>% filter(!is.na(ensembl_gene_id))
ratite_spec_genes_strict <- ratite_spec_genes %>% filter(entrezgene %in% ratite_spec_genes_gt$entrezgene) %>% filter(!is.na(ensembl_gene_id))

ratite_genes <- bsrel %>% filter(tsel.s >= 1, tsel.s > nsel.s, species_tree == "True")  %>% filter(!is.na(ensembl_gene_id))
ratite_genes_gt <- bsrel %>% filter(tsel.s >= 1, tsel.s > nsel.s, species_tree == "False")  %>% filter(!is.na(ensembl_gene_id))
ratite_genes_strict <- ratite_genes %>% filter(entrezgene %in% ratite_genes_gt$entrezgene)  %>% filter(!is.na(ensembl_gene_id))

background <- bsrel %>% distinct(entrezgene, .keep_all=TRUE)  %>% filter(!is.na(ensembl_gene_id))


bsrel_s_mf <- enrichGO(ratite_spec_genes_strict$ensembl_gene_id,'org.Gg.eg.db',pvalueCutoff = 0.05, universe=background$ensembl_gene_id,keyType="ENSEMBL",ont="MF") %>% as.data.frame
bsrel_s_bp <- enrichGO(ratite_spec_genes_strict$ensembl_gene_id,'org.Gg.eg.db',pvalueCutoff = 0.05, universe=background$ensembl_gene_id,keyType="ENSEMBL",ont="BP") %>% as.data.frame

bsrel_all_mf <- enrichGO(ratite_genes_strict$ensembl_gene_id,'org.Gg.eg.db',pvalueCutoff = 0.05, universe=background$ensembl_gene_id,keyType="ENSEMBL",ont="MF") %>% as.data.frame
bsrel_all_bp <- enrichGO(ratite_genes_strict$ensembl_gene_id,'org.Gg.eg.db',pvalueCutoff = 0.05, universe=background$ensembl_gene_id,keyType="ENSEMBL",ont="BP") %>% as.data.frame

bsrel_s_kegg <- enrichKEGG(ratite_spec_genes_strict$entrezgene, 'gga', pvalueCutoff = 0.05, universe = background$entrezgene, keyType = "ncbi-geneid")
bsrel_all_kegg <- enrichKEGG(ratite_genes_strict$entrezgene, 'gga', pvalueCutoff = 0.05, universe = background$entrezgene, keyType = "ncbi-geneid")

bsrel_s_mf
bsrel_s_bp
bsrel_all_mf
bsrel_all_bp
bsrel_s_kegg
bsrel_all_kegg


#work with tibbles
merged <- tbl_df(merged)
merged <- hog_info %>% dplyr::select(hog, dup_ct, missing_ct) %>% right_join(., merged) %>% tbl_df %>% 
  mutate(totbranch = tnon + nnon + tsel.n + nsel.n) %>% 
  mutate(total_sel.s = tsel.s + nsel.s, total_sel.n = tsel.n + nsel.n, target_prop.s = tsel.s/total_sel.s, target_prop.n = tsel.n / total_sel.n, target_lin = tsel.n + tnon) %>% 
  dplyr::select(-class)

## THIS CODE CHECKS FOR MISSING HOGS TO RERUN
#get missing and set up reruns
check_for_missing = subset(merged, select=c("hog", "tree", "treenum", "species_tree", "totbranch"))

table(check_for_missing$treenum, check_for_missing$tree, useNA="ifany")
check_for_missing=subset(check_for_missing, !is.na(treenum))
table(check_for_missing$totbranch==0, useNA='ifany')
#missing set to rerun
rerun<-subset(check_for_missing, totbranch==0 | is.na(totbranch))
#rerun rule: if tree1 is a rerun & the species tree, rerun tree1 and tree2
#if tree1 is a rerun and not the species tree, no tree2 to rerun
#if tree2 is a rerun and not the species tree, no need to rerun tree1
#however because of checks in script, should be able to just run straight up
write.table(unique(rerun$hog), col.names = F, row.names = F, quote=F, file="hogs_to_rerun_Nov2017")
#####

#ANALYSIS 
bsrel<-merged %>% filter(dup_ct == 0, species_tree == TRUE, totbranch > 0, missing_ct <= 5)
bsrel %>% distinct(hog) %>% summarize(count=n())

#test whether ratite-specific selection has functional enrichments
ratite_spec_genes <- bsrel %>% filter(tsel.s >= 1, nsel.s == 0) %>% dplyr::select(hog) %>% inner_join(., hog_to_gene) %>% dplyr::select(gene)
ratite_conv_genes <- bsrel %>% filter(tsel.s > 1) %>% dplyr::select(hog) %>% inner_join(., hog_to_gene) %>% dplyr::select(gene)
ratite_conv_spec_genes <- ratite_conv_genes %>% filter(gene %in% ratite_spec_genes$gene)


bsrel_sc_mf <- enrichGO(ratite_conv_spec_genes$gene,'org.Gg.eg.db',pvalueCutoff = 0.1, universe=backgroundset$gene,keytype="SYMBOL",ont="MF")
bsrel_sc_bp <- enrichGO(ratite_conv_spec_genes$gene,'org.Gg.eg.db',pvalueCutoff = 0.1, universe=backgroundset$gene,keytype="SYMBOL",ont="BP")

bsrel_c_mf <- enrichGO(ratite_conv_genes$gene,'org.Gg.eg.db',pvalueCutoff = 0.1, universe=backgroundset$gene,keytype="SYMBOL",ont="MF")
bsrel_c_bp <- enrichGO(ratite_conv_genes$gene,'org.Gg.eg.db',pvalueCutoff = 0.1, universe=backgroundset$gene,keytype="SYMBOL",ont="BP")

bsrel_s_mf <- enrichGO(ratite_spec_genes$gene,'org.Gg.eg.db',pvalueCutoff = 0.1, universe=backgroundset$gene,keytype="SYMBOL",ont="MF")
bsrel_s_bp <- enrichGO(ratite_spec_genes$gene,'org.Gg.eg.db',pvalueCutoff = 0.1, universe=backgroundset$gene,keytype="SYMBOL",ont="BP")

bsrel_sc_mf 
bsrel_sc_bp
bsrel_c_mf 
bsrel_c_bp
bsrel_s_mf
bsrel_s_bp

#convergence testing
barplot(table(bsrel$tsel.s[bsrel$nsel.s == 0 & bsrel$tsel.s >= 1]), ylab="Number of Genes Uniquely Selected in >1 Target Lineages", legend=T)
bsrel$prob_sel = (bsrel$total_sel.s / bsrel$totbranch)

compute_convergence <- function(prob = prob, target = target, total = total) {
  target<-rbinom(1, target, prob)
  nontarget<-rbinom(1, total-target, prob)
  if (nontarget == 0) {
    return(target)
  }
  else {
    return(NA_real_)
  }
}

exp_multiratite<-list()

for (i in 1:1000) {
  res<-mapply(compute_convergence, bsrel$prob_sel, bsrel$target_lin, bsrel$totbranch)
exp_multiratite[[i]]<-as.data.frame(table(res), stringsAsFactors = F)
}

perm_res <- bind_rows(exp_multiratite, .id="rep") %>% spread(res, Freq, fill=0, drop=FALSE) %>% 
  dplyr::rename(ct0 = `0`, ct1 = `1`, ct2 = `2`, ct3 = `3`, ct4=`4`, ct5=`5`) %>%
  mutate(selected = ct1+ct2+ct3+ct4+ct5, convergent = ct2+ct3+ct4+ct5, rel_conv = convergent/selected) %>% tbl_df

real_res <- as.data.frame(table(bsrel$tsel.s[bsrel$nsel.s==0])) %>% tbl_df %>% spread(Var1, Freq) %>%
  dplyr::rename(ct0 = `0`, ct1 = `1`, ct2 = `2`, ct3 = `3`) %>%
  mutate(ct4=0, ct5=0, selected = ct1+ct2+ct3+ct4+ct5, convergent = ct2+ct3+ct4+ct5, rel_conv = convergent/selected) %>% tbl_df

1-(sum(real_res$selected > perm_res$selected)/1000)
1-(sum(real_res$convergent > perm_res$convergent)/1000)
1-(sum(real_res$rel_conv > perm_res$rel_conv)/1000)

bsrel %>% filter(tsel.s > 2, nsel.s == 0) %>% dplyr::select(hog) %>% inner_join(., hog_to_gene) %>% arrange(gene) %>% print.data.frame

library(ggthemes)
#figure 2b
perm_res %>% ggplot(aes(x=rel_conv)) + theme_tufte() + geom_density() + labs(x="Proportion with convergent selection") +
  geom_vline(xintercept=real_res$rel_conv, col="red")
