setwd("~/Projects/birds/ratite_compgen/ratite-genomics/06_protein_coding_analysis/hyphy_bsrel/")
library(tidyverse)
library(clusterProfiler)
library(biomaRt)
library(ape)
library(gridExtra)

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
hog_final_key %>% filter(species_tree == "True") %>% count(hog) %>% filter(n==2)

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

count_targets <- function(nodelist, targets, return_full = FALSE) {
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
  
  if (return_full) {
    return(node_vec %>% filter(target_frac >= 1))
  }
  
  sum(node_vec$target_frac >= 1)
}

count_ratite_clades <- function(nodelist) {
  
  ratite_clades <- tibble(node = c("anoDid", "strCam", "rheAme", "rhePen", "rheAme-rhePen", "droNov", "casCas", "casCas-droNov", "aptHaa", "aptRow", "aptOwe", "aptHaa-aptOwe", "aptHaa-aptOwe-aptRow"),
                          clade = c("moa", "ostrich", "rhea", "rhea", "rhea", "emu", "emu", "emu", "kiwi", "kiwi", "kiwi", "kiwi", "kiwi"))
  
  ratite_counts <- count_targets(nodelist, ratite_clades$node, TRUE)
  
  inner_join(ratite_clades, ratite_counts, by=c("node" = "node")) %>% distinct(clade) %>% pull(clade) %>% paste0(collapse = "-")

}

ratite_tips = c("anoDid", "strCam", "rheAme", "rhePen", "droNov", "casCas", "aptHaa", "aptRow", "aptOwe")
corm_tip = c("nanHar")

bsrel_orig <- bsrel_orig %>% 
  mutate(ratite_selected_nom = map_int(nom_br, count_targets, ratite_tips)) %>%
  mutate(ratite_selected_str = map_int(strict_br, count_targets, ratite_tips)) %>%
  mutate(ratite_nom_clades = map_chr(nom_br, count_ratite_clades)) %>%
  mutate(ratite_str_clades = map_chr(strict_br, count_ratite_clades))

bsrel_ext <- bsrel_ext %>%
  mutate(ratite_selected_nom = map_int(nom_br, count_targets, ratite_tips)) %>%
  mutate(ratite_selected_str = map_int(strict_br, count_targets, ratite_tips)) %>%
  mutate(ratite_nom_clades = map_chr(nom_br, count_ratite_clades)) %>%
  mutate(ratite_str_clades = map_chr(strict_br, count_ratite_clades)) %>%
  mutate(corm_selected_nom = map_int(nom_br, count_targets, corm_tip)) %>%
  mutate(corm_selected_str = map_int(strict_br, count_targets, corm_tip))

bsrel_orig_clean <- bsrel_orig %>% filter(species_tree == "True", !is.na(total), total > 0, !grepl("XM", nom_br), !grepl("00", nom_br)) %>% distinct(hog, chr, entrezgene, .keep_all = TRUE)

bsrel_ext_clean <- bsrel_ext %>% filter(species_tree == "True", !is.na(total), total > 0, !grepl("XM", nom_br), !grepl("00", nom_br)) %>% distinct(hog, chr, entrezgene, .keep_all = TRUE)

pos_ext <- bsrel_ext_clean %>% dplyr::select(hog, entrezgene, sym = external_gene_name, strict_br) %>% 
  filter(!grepl("XM", strict_br))

ext_tree <- read.tree("../final_tree_ext_proteins.nwk")

all_tips <- ext_tree$tip.label[1:44]
all_tips_internal <- bsrel_ext_clean %>% distinct(nom_br) %>% pull(nom_br) %>% str_split(., coll(":")) %>% unlist %>% unique

count_selected_tips <- function(branch_string, tips, drop_internal = TRUE) {
  nodes <- str_split(branch_string, coll(":")) %>% unlist
  if (drop_internal) {
    nodes <- nodes[!grepl("-", nodes)]
  }
  selected<-as.numeric(tips %in% nodes)
  names(selected) <- tips
  return(selected)
}

pos_ext_wide <- apply(pos_ext[,4], 1, function(x) count_selected_tips(x, all_tips)) %>% t() %>% as.tibble() %>% 
  mutate(total_sel = rowSums(.), ratite_sel = anoDid+strCam+rheAme+rhePen+droNov+casCas+aptHaa+aptRow+aptOwe)

write_tsv(pos_ext_wide, path="bsrel_matrix.tsv")

##ANALYSIS##

bsrel_ext_clean %>% filter(ratite_selected_str == strict, strict > 0)
bsrel_ext_clean %>% filter(ratite_selected_str == strict, strict > 0) %>% count(ratite_str_clades) %>% print.data.frame
bsrel_ext_clean %>% filter(ratite_selected_str >= 1, ratite_selected_str/strict > 0.50) %>% count(ratite_str_clades) %>% print.data.frame

bsrel_ext_clean %>% filter(ratite_selected_str >= 1, ratite_selected_str/strict > 0.50, ratite_str_clades == "moa-ostrich-rhea-emu-kiwi")

background_orig <- bsrel_orig_clean %>% distinct(entrezgene)
background_ext <- bsrel_ext_clean %>% distinct(entrezgene)

##CONVERGENCE##

#ratite vs cormorant

pos_ext_wide %>% colSums(.)
pos_ext_wide %>% filter(total_sel == ratite_sel + 1, ratite_sel > 0) %>% colSums(.) / colSums(pos_ext_wide)
bsrel_ext_clean %>% filter(ratite_selected_str >= 1, ratite_selected_str + 1 == strict, corm_selected_str == 1)

bsrel_ext_clean <- bsrel_ext_clean %>% mutate(prob_sel_str = strict/total, prob_sel_nom = nom/total)
bsrel_orig_clean <- bsrel_orig_clean %>% mutate(prob_sel_str = strict/total, prob_sel_nom = nom/total)

#run this once then load data

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

conv_perms_orig<-list()
conv_perms_ext<-list()
conv_perms_trio<-list()

for (i in 1:10000) {
  res<-mapply(compute_convergence, bsrel_orig_clean$prob_sel_str, 4, bsrel_orig_clean$total-4)
  res2<-mapply(compute_convergence, bsrel_ext_clean$prob_sel_str, 5, bsrel_ext_clean$total-4)
  res3<-mapply(compute_convergence, bsrel_ext_clean$prob_sel_str, 3, bsrel_ext_clean$total-6)
  conv_perms_orig[[i]]<-as.data.frame(table(res), stringsAsFactors = F)
  conv_perms_ext[[i]]<-as.data.frame(table(res2), stringsAsFactors = F)
  conv_perms_trio[[i]]<-as.data.frame(table(res3), stringsAsFactors = F)
}

perm_res_orig <- bind_rows(conv_perms_orig, .id="rep") %>% spread(res, Freq, fill=0, drop=FALSE) %>% 
  dplyr::rename(ct0 = `0`, ct1 = `1`, ct2 = `2`, ct3 = `3`, ct4 = `4`) %>%
  mutate(selected = ct1+ct2+ct3+ct4, convergent = ct2+ct3+ct4, rel_conv = convergent/selected) %>% tbl_df

perm_res_ext <- bind_rows(conv_perms_ext, .id="rep") %>% spread(res2, Freq, fill=0, drop=FALSE) %>% 
  dplyr::rename(ct0 = `0`, ct1 = `1`, ct2 = `2`, ct3 = `3`, ct4 = `4`, ct5 = `5`) %>%
  mutate(selected = ct1+ct2+ct3+ct4+ct5, convergent = ct2+ct3+ct4+ct5, rel_conv = convergent/selected) %>% tbl_df

perm_res_trio <- bind_rows(conv_perms_trio, .id="rep") %>% spread(res3, Freq, fill=0, drop=FALSE) %>% 
  dplyr::rename(ct0 = `0`, ct1 = `1`, ct2 = `2`, ct3 = `3`) %>%
  mutate(selected = ct1+ct2+ct3, convergent = ct2+ct3, rel_conv = convergent/selected) %>% tbl_df

##end run once

write_tsv(perm_res_trio, path="perm_res_trio.tsv")
write_tsv(perm_res_ext, path="perm_res_ext.tsv")
write_tsv(perm_res_orig, path="perm_res_orig.tsv")

perm_res_ext <- read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/06_protein_coding_analysis/hyphy_bsrel/perm_res_ext.tsv")
perm_res_trio <- read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/06_protein_coding_analysis/hyphy_bsrel/perm_res_trio.tsv")
perm_res_orig <- read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/06_protein_coding_analysis/hyphy_bsrel/perm_res_orig.tsv")

#real results

orig_counts <- bsrel_orig_clean %>% 
  filter(ratite_selected_str >= 1,ratite_selected_str == strict, species_tree == "True", ratite_str_clades != "") %>% 
  count(ratite_str_clades) %>% mutate(clade_ct = str_count(ratite_str_clades, coll("-"))+1)

ext_counts <- bsrel_ext_clean %>% 
  filter(ratite_selected_str >= 1,ratite_selected_str == strict, species_tree == "True", ratite_str_clades != "") %>% 
  count(ratite_str_clades) %>% mutate(clade_ct = str_count(ratite_str_clades, coll("-"))+1)

orig_counts %>% group_by(clade_ct) %>% summarize(sum(n))
ext_counts %>% group_by(clade_ct) %>% summarize(sum(n))

#extended counts P-value:

(1+sum(perm_res_ext$rel_conv >= (50+11+2)/(50+11+2+175)))/10000
(1+sum(perm_res_orig$rel_conv >= (1+26)/(1+26+325)))/10000

#trios

#convergent = all three lineages, so moa-ostrich and one more
#get mean of convergence with each third lineage:
ext_counts %>% filter(grepl("kiwi", ratite_str_clades))
kiwi<-(2+6+2+1+1+1)/(ext_counts %>% filter(grepl("kiwi", ratite_str_clades)) %>% pull(n) %>% sum(.))
ext_counts %>% filter(grepl("rhea", ratite_str_clades))
rhea<-(1+3+2+1+4+1+1)/(ext_counts %>% filter(grepl("rhea", ratite_str_clades)) %>% pull(n) %>% sum(.))
ext_counts %>% filter(grepl("emu", ratite_str_clades))
emu<-(14+2+3+2+2+1)/(ext_counts %>% filter(grepl("emu", ratite_str_clades)) %>% pull(n) %>% sum(.))

#plotting

plot_orig <- perm_res_orig %>% ggplot(aes(rel_conv)) + geom_histogram(bins=50, fill="steelblue") + theme_minimal(base_size = 14) + xlab("Proportion Selected Genes with Evidence for Convergence") + coord_cartesian(xlim=c(0,0.30)) + geom_vline(xintercept=0.07670455, col="red")
plot_ext <- perm_res_ext %>% ggplot(aes(rel_conv)) + geom_histogram(bins=50, fill="steelblue") + theme_minimal(base_size = 14) + xlab("Proportion Selected Genes with Evidence for Convergence") + coord_cartesian(xlim=c(0,0.30)) + geom_vline(xintercept=0.2647059, col="red")

figS9<-grid.arrange(plot_ext, plot_orig)
ggsave("~/Projects/birds/ratite_compgen/manuscript/ScienceSubmissionRev1/FullDraftDec11/FigS9-Dec11.pdf", figS9)

#checking consistency of results across original and extended runs:
merged <- full_join(bsrel_ext_clean, bsrel_orig_clean, by=c("hog" = "hog"), suffix=c(".ext", ".org"))

merged %>% ggplot(aes(nom.ext, nom.org)) + geom_point(alpha=0.1)
merged %>% ggplot(aes(strict.org, nom.org)) + geom_point(alpha=0.1)

merged %>% filter(strict.ext > 40, strict.org < 2) %>% View()

#enrichments for convergent genes? what are the ratite-specific and convergent genes?

ratite_selected_genes <- bsrel_ext_clean %>% 
  filter(ratite_selected_str >= 1,ratite_selected_str == strict, species_tree == "True", ratite_str_clades != "") %>% 
  mutate(clade_ct = str_count(ratite_str_clades, coll("-"))+1) %>% 
  dplyr::select(hog, chr, entrezgene, external_gene_name, strict, strict_br, ratite_selected_str, ratite_str_clades, clade_ct)

#GO / KEGG#

ratite_spec_genes_orig <- bsrel_orig_clean %>% filter(ratite_selected_str >= 1,ratite_selected_str == strict, species_tree == "True")
ratite_spec_genes_ext <- bsrel_ext_clean %>% filter(ratite_selected_str >= 1,ratite_selected_str == strict, species_tree == "True")

ratite_enrich_genes_ext <- bsrel_ext_clean %>% filter(ratite_selected_str >= 1,ratite_selected_str/strict > 0.50, species_tree == "True")

ratite_enrich_genes_orig <- bsrel_orig_clean %>% filter(ratite_selected_str >= 1,ratite_selected_str/strict > 0.50, species_tree == "True")

background_ext <- bsrel_ext_clean %>% filter(species_tree == "True")
background_orig <- bsrel_orig_clean %>% filter(species_tree == "True")

#Ratite-specific
enrichGO(ratite_spec_genes_ext$entrezgene,'org.Gg.eg.db',pvalueCutoff = 0.20, universe=background_ext$entrezgene,keyType="ENTREZID",ont="MF") %>% as.data.frame %>% dplyr::select(ID:qvalue)
enrichGO(ratite_spec_genes_ext$entrezgene,'org.Gg.eg.db',pvalueCutoff = 0.20, universe=background_ext$entrezgene,keyType="ENTREZID",ont="BP") %>% as.data.frame %>% dplyr::select(ID:qvalue)
enrichKEGG(ratite_spec_genes_ext$entrezgene,'gga',pvalueCutoff = 0.20, universe=background_ext$entrezgene,keyType="ncbi-geneid") %>% as.data.frame %>% dplyr::select(ID:qvalue)

#Ratite enriched
enrichGO(ratite_enrich_genes_ext$entrezgene,'org.Gg.eg.db',pvalueCutoff = 0.20, universe=background_ext$entrezgene,keyType="ENTREZID",ont="MF") %>% as.data.frame %>% dplyr::select(ID:qvalue)
summary(enrichGO(ratite_enrich_genes_ext$entrezgene,'org.Gg.eg.db',pvalueCutoff = 0.20, universe=background_ext$entrezgene,keyType="ENTREZID",ont="BP")) %>% as.data.frame %>% dplyr::select(ID:qvalue) 
enrichKEGG(ratite_enrich_genes_ext$entrezgene,'gga',pvalueCutoff = 0.20, universe=background_ext$entrezgene,keyType="ncbi-geneid") %>% as.data.frame %>% dplyr::select(ID:qvalue)

#Ratite convergent
enrichGO(ratite_selected_genes$entrezgene[ratite_selected_genes$clade_ct > 1],'org.Gg.eg.db',pvalueCutoff = 0.20, universe=background_ext$entrezgene,keyType="ENTREZID",ont="BP")
enrichGO(ratite_selected_genes$entrezgene[ratite_selected_genes$clade_ct > 1],'org.Gg.eg.db',pvalueCutoff = 0.20, universe=background_ext$entrezgene,keyType="ENTREZID",ont="MF")
enrichKEGG(ratite_selected_genes$entrezgene[ratite_selected_genes$clade_ct > 1],'gga',pvalueCutoff = 0.20, universe=background_ext$entrezgene,keyType="ncbi-geneid")



