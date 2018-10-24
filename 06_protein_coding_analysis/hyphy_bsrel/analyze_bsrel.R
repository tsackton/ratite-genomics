setwd("~/Projects/birds/ratite_compgen/ratite-genomics/06_protein_coding_analysis/hyphy_bsrel/")
library(tidyverse)
library(clusterProfiler)
library(biomaRt)

##ANNOTATION -- code copied from https://github.com/ajshultz/avian-immunity/blob/master/ComparativeGenomics/PAML_Analyze/02_GeneID_Annotation.R ##

#Get CDS <-> GENE key for chicken
gg4_ft<-read_tsv("../../04_wga/03_ce_annotation/GCF_000002315.3_Gallus_gallus-4.0_genes_table.txt", col_types = "ccccciicccc") 
names(gg4_ft) = c("feature", "class", "asm", "chr", "genomic_acc", "start", "end", "accession", "nr_ref", "sym", "entrezgene")

#Get ACC from actual PAML alignments
hog_to_acc <- read_tsv("../../03_homology/HOG_final_alignment_seqids") %>% 
  filter(Taxon == "galGal") %>%
  mutate(hog = sub("HOG2_", "", HOG)) %>%
  dplyr::select(hog, trans_acc = Transcript)

#merge with feature table
gg4_ft <- gg4_ft %>% 
  filter(!is.na(accession)) %>% 
  right_join(hog_to_acc, by=c("accession" = "trans_acc")) %>%
  dplyr::select(hog, chr, entrezgene, sym)

####Add human NCBI gene IDs from biomaRt
ensembl = useMart("ensembl")

#Use correct datasets - sometimes throws an error for some reason saying the given dataset is not valid. However, it works if you try again once or twice
ensembl_gg = useDataset("ggallus_gene_ensembl",mart=ensembl)
ensembl_hs = useDataset("hsapiens_gene_ensembl",mart=ensembl)

#Get human ensembl IDs for homologs
human_ids_gg <- getBM(attributes=c("ensembl_gene_id","external_gene_name","hsapiens_homolog_ensembl_gene","hsapiens_homolog_orthology_confidence"),
                      filters=c("with_hsapiens_homolog"),
                      values=TRUE,
                      mart=ensembl_gg) %>% tbl_df

#Get ensembl ID to NCBI entrezgene mapping for three different species
gg_ens_to_entrezgene <- getBM(attributes=c("ensembl_gene_id","entrezgene"),
                              mart=ensembl_gg) %>%
  tbl_df %>%
  filter(!is.na(entrezgene)) %>%
  mutate(entrezgene = as.character(entrezgene))

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

#make final hog key
hog_key_gg <- gg4_ft %>% left_join(gg_trans_table, by=c("entrezgene" = "entrezgene"))

full_list <- read_tsv("../paml_ancrec/paml_M0_parsed.txt.gz", col_types = "cccc?????") %>% 
  filter(!is.na(hog), model == "ancrec") %>%
  mutate(tree = paste0("tree", treenum)) %>%
  dplyr::select(hog, tree, species_tree)

#make a simpler hog key with just chicken entrezgene ensembl_gene and external gene name
hog_final_key <- hog_key_gg %>% full_join(full_list, by=c("hog" = "hog")) %>%
  distinct(hog, tree, species_tree, chr, entrezgene, ensembl_gene_id, external_gene_name) %>% 
  filter(!is.na(tree))

hog_final_key %>% filter(species_tree == "False") %>% count(hog) %>% count(n)
hog_final_key %>% filter(species_tree == "True") %>% count(hog) %>% count(n)

hog_final_key %>% count(hog) %>% filter(n > 2)
hog_final_key %>% filter(species_tree == "True") %>% count(hog) %>% filter(n >= 2)
hog_final_key %>% filter(hog == 27140)

#hog final key still has some duplication - looks like mostly cases where an entrez gene id maps to two ENS gene ids

##LOAD DATA## 
bsrel_cols <- c("tree", "hog", "strict", "nom", "total", "strict_br", "nom_br", "newick")
bsrel_ext <- read_tsv("bsrel_res_parsed_extended_2018-10-22.txt", col_names = bsrel_cols)
bsrel_orig <- read_tsv("bsrel_res_parsed_original_2018-10-17.txt", col_names = bsrel_cols)

##ADD GENES, CLEAN UP##

bsrel_ext <- hog_final_key %>% distinct(hog, tree, entrezgene, .keep_all = TRUE) %>%
  mutate(hog = as.integer(hog)) %>%
  right_join(bsrel_ext, by=c("hog" = "hog", "tree" = "tree")) %>%
  filter(!is.na(species_tree))

bsrel_orig <- hog_final_key %>% distinct(hog, tree, entrezgene, .keep_all = TRUE) %>%
  mutate(hog = as.integer(hog)) %>%
  right_join(bsrel_orig, by=c("hog" = "hog", "tree" = "tree")) %>%
  filter(!is.na(species_tree))


#this looks very good now, ready for analysis

#first need to write function to take tree string and compute # targets vs # non-targets

count_targets <- function(nodelist, targets) {
  #nodelist will need to be split on ":"
  #then each element needs to be checked for presence of targets
  #and internal nodes computed to assign to target or non-target (tips are easy)
  
  ##INTERNAL FUNCS##
  #parse a single node (can be internal) into pieces, fraction of pieces that are in target
  compute_target_frac <- function(node, targets) {
    test<-strsplit(node, "-", fixed=TRUE) %>% unlist
    target_hits <- sum(test %in% targets)
    total_len <- length(test)
    return (target_hits/total_len)
  }
  
  #first split nodelist
  node_vec <- tibble(node = strsplit(nodelist, ":", fixed=TRUE) %>% unlist) %>% 
    group_by(node) %>% 
    mutate(target_frac = compute_target_frac(node, targets))
  
  sum(node_vec$target_frac >= 1)
}

ratite_tips = c("anoDid", "strCam", "rheAme", "rhePen", "droNov", "casCas", "aptHaa", "aptRow", "aptOwe")
corm_tip = c("nanHar")

bsrel_orig <- bsrel_orig %>% 
  mutate(ratite_selected_nom = map_int(nom_br, count_targets, ratite_tips)) %>%
  mutate(ratite_selected_str = map_int(strict_br, count_targets, ratite_tips))

bsrel_ext <- bsrel_ext %>%
  mutate(ratite_selected_nom = map_int(nom_br, count_targets, ratite_tips)) %>%
  mutate(ratite_selected_str = map_int(strict_br, count_targets, ratite_tips)) %>%
  mutate(corm_selected_nom = map_int(nom_br, count_targets, corm_tip)) %>%
  mutate(corm_selected_str = map_int(strict_br, count_targets, corm_tip))

##ANALYSIS##

ratite_spec_genes_orig <- bsrel_orig %>% filter(ratite_selected_str >= 1,ratite_selected_str == strict, species_tree == "True")
ratite_spec_genes_ext <- bsrel_ext %>% filter(ratite_selected_str >= 1,ratite_selected_str == strict, species_tree == "True")



### OLD CODE BELOW ###


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
