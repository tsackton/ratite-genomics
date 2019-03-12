setwd("~/Projects/birds/ratite_compgen/ratite-genomics/06_protein_coding_analysis/hyphy_relax/")
library(data.table)

#ratites
rp<-read.table("relax_pval.all", header=F)
rk<-read.table("relax_K.all", header=F, fill=T)

rp=rp[,c(1,4,5,6,8)]
rk=rk[,c(1,4,5,6,8)]

names(rp)=c("set", "hog", "tree", "type", "pval")
names(rk)=c("set", "hog", "tree", "type", "K")
relax=merge(rp,rk)

relax=subset(relax, pval != "-------" & !is.na(K))
relax$pval=as.numeric(as.character(relax$pval))
relax=unique(relax)

#finished, good runs in relax now

#check for missing
ancrec.parsed<-fread("gunzip -c ../paml_ancrec/paml_M0_parsed.txt.gz")
ancrec.treekey<-ancrec.parsed[,c("hog", "treenum", "species_tree"), with=FALSE]
hog_info<-read.table("../../03_homology/all_hog_info.tsv", sep="\t", header=T)

hog_info$has_species_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree]
hog_info$has_gene_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree == F]

#merge tree info with relax
relax$tree = as.integer(sub("tree", "", relax$tree))
relax = merge(relax, ancrec.treekey, by.x=c("hog", "tree"), by.y=c("hog", "treenum"))
relax.sp = subset(relax, species_tree)
relax.gt = subset(relax, species_tree==FALSE)

#convert to data table
relax = as.data.table(relax)

hog_counts<-relax[,.N,by=.(hog,set)]
hogs_to_run = hog_info$hog[hog_info$has_species_tree]

write.table(hogs_to_run[!(hogs_to_run %in% hog_counts$hog[hog_counts$set=="Ratite" & hog_counts$N == 2])], file="ratite_reruns_Dec2017", quote=F, row.names = F, col.names = F)
write.table(hogs_to_run[!(hogs_to_run %in% hog_counts$hog[hog_counts$set=="VL" & hog_counts$N == 2])], file="vl_reruns_Dec2017", quote=F, row.names = F, col.names = F)
write.table(hogs_to_run[!(hogs_to_run %in% hog_counts$hog[hog_counts$set=="RND" & hog_counts$N == 2])], file="rand_reruns_Dec2017", quote=F, row.names = F, col.names = F)

## ANALYSIS

#make subset
hog_counts_all <- relax[,.N,by=hog]
hogs_to_use = hog_counts_all$hog

relax<-merge(relax, hog_info, by="hog")
missing_cutoff = 2
relax.clean = subset(relax, (hog %in% hogs_to_use) & dup_ct == 0 & missing_ct <= missing_cutoff, select=c("hog", "set", "type", "pval", "K"))

#add qvalues
relax.clean[,qval := p.adjust(pval, method="fdr"), by=.(set, type)]
relax.clean$sample = paste0(relax.clean$set, ".", relax.clean$type)
table(relax.clean$qval < 0.05, relax.clean$sample)

#add sig key
relax.clean$sig = as.numeric(relax.clean$qval < 0.05) * sign(1 - relax.clean$K)
#positive is increased selection, negative is relaxed selection

table(relax.clean$sig, relax.clean$sample)

#compare with other branch tests

vl_rer <- res_RER_all %>% filter(targets == "vl", tree == "red") %>% mutate(hog = as.integer(sub("HOG", "", gene)))
rat_rer <- res_RER_all %>% filter(targets == "rat", tree == "red") %>% mutate(hog = as.integer(sub("HOG", "", gene)))
vl.cmp = relax.clean %>% filter(sample == "VL.all") %>% inner_join(vl_rer, by=c("hog"="hog"), suffix=c(".relax", ".btest"))
ratite.cmp = relax.clean %>% filter(sample == "Ratite.all") %>% inner_join(rat_rer, by=c("hog"="hog"), suffix=c(".relax", ".btest"))

#extended figure 6
vl.cmp %>% ggplot(aes(x=-log10(pval), y=-log10(P))) +
  geom_hex(bins=50) + geom_hex(bins=50) + labs(x="HyPhy RELAX (rank transformed)", y="Chikina et al (rank transformed)")
cor.test((1-vl.cmp$K), vl.cmp$est, method="sp")

ratite.cmp %>% mutate(hy_rank = rank(1-K), bt_rank = rank(est)) %>% ggplot(aes(x=hy_rank, y=bt_rank)) +
  geom_hex(bins=50) + geom_hex(bins=50) + labs(x="HyPhy RELAX (rank transformed)", y="Chikina et al (rank transformed)")

cor.test((1-ratite.cmp$K), ratite.cmp$est, method="sp")

#extended table 6
table(vl.cmp$p.adj < 0.05 & vl.cmp$direction == "up", vl.cmp$sig == -1) %>% fisher.test
table(vl.cmp$p.adj < 0.05 & vl.cmp$direction == "down", vl.cmp$sig == 1) %>% fisher.test

table(ratite.cmp$p.adj < 0.05 & ratite.cmp$direction == "up", ratite.cmp$sig == -1) %>% fisher.test
table(ratite.cmp$p.adj < 0.05 & ratite.cmp$direction == "down", ratite.cmp$sig == 1) %>% fisher.test


###FINAL PLOTS##


ratite.cmp %>% filter(abs(log10(K)) <= 1) %>% ggplot(aes(Rho, log10(K))) +  geom_point(alpha=0.25) + labs(x="Effect size (Rho) from RERconverge", y="Effect size (log K) from RELAX") + theme_classic(base_size = 14) + geom_smooth(method="lm", se=FALSE, col="red")
ggsave("~/Projects/birds/ratite_compgen/manuscript/ScienceSubmissionRev1/FullDraftDec11/FigS7-Dec11.pdf")
cor.test(ratite.cmp$K, ratite.cmp$Rho, method="sp")

