#branch permutation processing

setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/paml_branch/")
library(tidyverse)
library(qvalue)
library(clusterProfiler)
library(org.Gg.eg.db)
library(ggthemes)

#load data -- random species permutations

rand<-read_tsv("~/Projects/birds/ratite_compgen/data/branch_perms/all_rand.merged2.gz", col_names=c("est", "pval", "hog", "tips", "perm"))
perm.counts.rand <- rand %>% group_by(perm) %>% summarize(tip=unique(tips), sig.up = sum(pval < 0.01 & est > 0, na.rm=T), sig.down = sum(pval < 0.01 & est < 0, na.rm=T)) %>% distinct(tip, .keep_all=TRUE)
rand %>% distinct(perm) %>% count()

ratite<-read_tsv("~/Projects/birds/ratite_compgen/data/branch_perms/ratite_rand.merged2.gz", col_names=c("est", "pval", "hog", "tips", "perm"))
perm.counts.ratite <- ratite %>% group_by(perm) %>% summarize(tip=unique(tips), sig.up = sum(pval < 0.01 & est > 0, na.rm=T), sig.down = sum(pval < 0.01 & est < 0, na.rm=T)) %>% distinct(tip, .keep_all=TRUE)
ratite %>% distinct(perm) %>% count()

vl<-read_tsv("~/Projects/birds/ratite_compgen/data/branch_perms/vl_rand.merged2.gz", col_names=c("est", "pval", "hog", "tips", "perm"))
perm.counts.vl <- vl %>% group_by(perm) %>% summarize(tip=unique(tips), sig.up = sum(pval < 0.01 & est > 0, na.rm=T), sig.down = sum(pval < 0.01 & est < 0, na.rm=T)) %>% distinct(tip, .keep_all=TRUE)  %>% filter(grepl("calAnn", tip), grepl("melUnd", tip))
vl %>% distinct(perm) %>% count()

#random trios analysis 
#starting again from input for figure

rand_counts <- perm.counts.rand %>% filter(!(tip %in% perm.counts.vl$tip), !(tip %in% perm.counts.ratite$tip)) 

(sum(rand_counts$sig.down >= mean(perm.counts.vl$sig.down))+1)/(length(rand_counts$sig.down)+1)
(sum(rand_counts$sig.down >= median(perm.counts.vl$sig.down))+1)/(length(rand_counts$sig.down)+1)

(sum(rand_counts$sig.up >= mean(perm.counts.vl$sig.up))+1)/(length(rand_counts$sig.up)+1)
(sum(rand_counts$sig.up >= median(perm.counts.vl$sig.up))+1)/(length(rand_counts$sig.up)+1)

(sum(rand_counts$sig.down >= mean(perm.counts.ratite$sig.down))+1)/(length(rand_counts$sig.down)+1)
(sum(rand_counts$sig.down >= median(perm.counts.ratite$sig.down))+1)/(length(rand_counts$sig.down)+1)

(sum(rand_counts$sig.up >= mean(perm.counts.ratite$sig.up))+1)/(length(rand_counts$sig.up)+1)
(sum(rand_counts$sig.up >= median(perm.counts.ratite$sig.up))+1)/(length(rand_counts$sig.up)+1)

#make plotting datasets
rand_plot <- rand %>% 
  filter(!(tips %in% perm.counts.vl$tip), !(tips %in% perm.counts.ratite$tip)) %>% 
  distinct(tips, pval, .keep_all = TRUE) %>% dplyr::select(pval, tips)  %>% 
  mutate(logp = log10(pval)) 

vl_plot <- vl %>% 
  filter(grepl("calAnn", tips), grepl("melUnd", tips)) %>%
  distinct(tips, pval, .keep_all = TRUE) %>% 
  dplyr::select(pval, tips) %>% 
  mutate(logp = log10(pval))
  
ratite_plot <- ratite %>% 
  distinct(tips, pval, .keep_all = TRUE) %>% 
  dplyr::select(pval, tips) %>% 
  mutate(logp = log10(pval))

p_ecdf <- rand_plot%>%
  ggplot(aes(x=logp, group=tips)) + 
  stat_ecdf(geom="step", pad=FALSE, na.rm=TRUE, col="gray50", alpha=0.1) + 
  stat_ecdf(data=vl_plot, geom="step", pad=FALSE, na.rm=TRUE, col="purple") +
  stat_ecdf(data=ratite_plot, geom="step", pad=FALSE, na.rm=TRUE, col="red")

p_ecdf
ggsave("ecdf_option1.pdf")

p_ecdf2 <- rand_plot %>%
  ggplot(aes(x=pval, group=tips)) + 
  geom_line(stat="density", col="gray50", alpha=1/50) + 
  geom_line(data=vl_plot, stat="density", col="purple") +
  geom_line(data=ratite_plot, stat="density", col="red")

p_ecdf2
ggsave("density_option2.pdf")

p_ecdf3 <- rand_plot %>%
  ggplot(aes(x=pval, group=tips)) + 
  geom_line(stat="ecdf", col="gray50", alpha=1/50) + 
  geom_line(data=vl_plot, stat="ecdf", col="purple") +
  geom_line(data=ratite_plot, stat="ecdf", col="red")

p_ecdf3
ggsave("ecdf_option3.pdf")

#FINAL PLOT
rand_plot %>%
  ggplot(aes(x=pval, group=tips)) + theme_classic() +
  geom_line(stat="ecdf", col="gray50", alpha=1/30) + 
  geom_line(data=vl_plot, stat="ecdf", col="purple") +
  geom_line(data=ratite_plot, stat="ecdf", col="red") +
  geom_vline(xintercept = 0.01)
ggsave("Figure2A_new.pdf")

#inset dataset
rand_counts %>% gather(key = "dir", value = "count", sig.up:sig.down) %>% 
  mutate(count = count * sign(0.1-as.numeric(dir=="sig.down"))) %>%
  ggplot(aes(x=count)) + geom_histogram(binwidth = 5, boundary=0, closed="left", fill="gray70") +
  theme_classic() + labs(x="Number of genes significant at nominal P < 0.01") + 
  geom_segment(aes(x=mean(perm.counts.ratite$sig.up), xend=mean(perm.counts.ratite$sig.up), y=50, yend=0), arrow=arrow(length=unit(0.1, "cm")), colour="red") + 
  geom_segment(aes(x=-1*mean(perm.counts.ratite$sig.down), xend=-1*mean(perm.counts.ratite$sig.down), y=50, yend=0), arrow=arrow(length=unit(0.1, "cm")), colour="red") +
  geom_segment(aes(x=mean(perm.counts.vl$sig.up), xend=mean(perm.counts.vl$sig.up), y=50, yend=0), arrow=arrow(length=unit(0.1, "cm")), colour="purple") + 
  geom_segment(aes(x=-1*mean(perm.counts.vl$sig.down), xend=-1*mean(perm.counts.vl$sig.down), y=50, yend=0), arrow=arrow(length=unit(0.1, "cm")), colour="purple") 
ggsave("Figure2A_inset.pdf")


#simple perms -- test for individuals hogs
vl_perm <- read_tsv("~/Projects/birds/ratite_compgen/data/branch_perms/simple_vl.merged.gz", col_names=c("est", "pval", "hog", "perm"))
ratite_perm <- read_tsv("~/Projects/birds/ratite_compgen/data/branch_perms/simple_ratite.merged.gz", col_names=c("est", "pval", "hog", "perm"))

#count number of successful permutatios
vl_perm %>% distinct(perm) %>% count()
ratite_perm %>% distinct(perm) %>% count()

#read real results
vl_real <- read_tsv("vl.real")
ratite_real <- read_tsv("ratite.real")

#add empirical p-values / fdr to real results

vl_results <- inner_join(vl_perm, vl_real, by=c("hog" = "hog"), suffix=c(".perm", ".real"))
vl_epval <- vl_results %>% group_by(hog) %>% summarize(est = unique(est.real), p.down = (sum(est.perm <= est.real, na.rm=T)+1)/(n()+1), p.up = (sum(est.perm >= est.real, na.rm=T)+1)/(n()+1))

ratite_results <- inner_join(ratite_perm, ratite_real, by=c("hog" = "hog"), suffix=c(".perm", ".real"))
ratite_epval <- ratite_results %>% group_by(hog) %>% summarize(est = unique(est.real), p.down = (sum(est.perm <= est.real, na.rm=T)+1)/(n()+1), p.up = (sum(est.perm >= est.real, na.rm=T)+1)/(n()+1))

#ratite_epval calcs
summary(qvalue(ratite_epval$p.down))
summary(qvalue(ratite_epval$p.up))
summary(qvalue(ratite_real$pval))

#vl_epval calcs
summary(qvalue(vl_epval$p.down))
summary(qvalue(vl_epval$p.up))
summary(qvalue(vl_real$pval))

#vl functional enrichment
#need to load all the necessary files for hog -> chicken results
#load hog <-> chicken info
hog_to_gene <- read.table("../HOG_final_alignment_seqids", header=T, stringsAsFactors =F)
hog_to_gene <- hog_to_gene %>% tbl_df %>%
  mutate(hog = as.integer(sub("HOG2_", "", HOG, fixed=T))) %>% 
  filter(Taxon == "galGal") %>%
  dplyr::select(-HOG, -Fasta_seqid, -Taxon, -Transcript) %>%
  group_by(hog) %>%
  mutate(gene = paste0(Gene, sep="", collapse=";")) %>%
  dplyr::select(-Gene)

vl_epval <- inner_join(vl_epval, hog_to_gene, by=c("hog"="hog"))
targetset.up <- vl_epval %>% mutate(qvalue = qvalue(p.up)$qvalues) %>% filter(qvalue <= 0.1) %>% dplyr::select(gene)
targetset.down <- vl_epval %>% mutate(qvalue = qvalue(p.down)$qvalues) %>% filter(qvalue <= 0.1) %>% dplyr::select(gene)
backgroundset <- vl_epval %>% dplyr::select(gene)

vl_mf.up <- enrichGO(targetset.up$gene,'org.Gg.eg.db',pvalueCutoff = 0.10, universe=backgroundset$gene,keytype="SYMBOL",ont="MF")
vl_mf.down <- enrichGO(targetset.down$gene,'org.Gg.eg.db',pvalueCutoff = 0.10, universe=backgroundset$gene,keytype="SYMBOL",ont="MF")
vl_bp.up <- enrichGO(targetset.up$gene,'org.Gg.eg.db',pvalueCutoff = 0.10, universe=backgroundset$gene,keytype="SYMBOL",ont="BP")
vl_bp.down <- enrichGO(targetset.down$gene,'org.Gg.eg.db',pvalueCutoff = 0.10, universe=backgroundset$gene,keytype="SYMBOL",ont="BP")

vl_mf.up
vl_mf.down
vl_bp.up
vl_bp.down



