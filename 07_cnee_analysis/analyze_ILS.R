# read gene tree result
#setwd("~/Phd_2/Bird")
#gene_tree <- read.table("gene_trees/palaeognath_CNEE_genetree_summary.txt", sep="\t",stringsAsFactors = F, header=T)
#rownames(gene_tree) <- gene_tree$Locus

#score <- read.table("Gain_1012/Combined_elem_lik.txt", header=T)
#Z <- read.table("Gain_1012/Combined_post_Z_M2.txt", header=T)

#target_species <- c("strCam","rhePen","rheAme","casCas","droNov","aptRow","aptHaa","aptOwe","anoDid")
#outgroup <- c("allMis","allSin","croPor","gavGan","chrPic","cheMyd","anoCar")
#target_col <- which(colnames(Z) %in% sapply(target_species, function(x) paste0(x,"_2")))
#out_col <- which(colnames(Z) %in% sapply(outgroup, function(x) paste0(x,"_2"))) #37-43
#conserve_col <- grep("_2", colnames(Z))[1:43]
#conserve_col <- setdiff(conserve_col, c(target_col, out_col))

#target_col <- which(colnames(Z) %in% sapply(target_species, function(x) paste0(x,"_3")))
#selZ <- apply(Z,1, function(x) all(x[conserve_col] >0.8) & any(x[target_col] >0.8))
##selZ0 <- apply(Z,1, function(x) all(x[1:36] == 1)) # all conserved


#score$AIC = gene_tree[score$ID, "deltaAIC"]
#score$RF = gene_tree[score$ID, "normalized_RF"]


#sel <- which(score$logBF1>=10 & score$logBF2 >= 1) 
#sel0 <-  which(score$logBF1>=10 & score$logBF2 >= 1) 

#selZ <- intersect(which(selZ), sel0)


#score$sel = 0
#score[sel, "sel"] = 1

#score$selZ = 0
#score[selZ, "selZ"] = 1

score <- read_tsv("05_ILS/score_gene_trees.txt")

wilcox.test(score$AIC ~ score$sel)
wilcox.test(score$AIC ~ score$selZ)
wilcox.test(score$RF ~ score$sel)
wilcox.test(score$RF ~ score$selZ)

#base graphics version

pdf("gene_trees/RAR_tree_boxplot1.pdf")
par(mfrow=c(2,2))
boxplot(AIC~selZ, score, notch=T, names = c("281362\nnon RAR", "2639 RAR\nBF1>=10, BF2>=1"), outline=F, ylab="AIC",col=c("deepskyblue2", "brown"))
boxplot(RF~sel, score, notch=T, names = c("281362\nnon RAR", "2639 RAR\nBF1>=10, BF2>=1"), outline=F, ylab="normalized RF",col=c("deepskyblue2", "brown"))
boxplot(AIC~selZ, score, notch=T, names = c("282709\nnon RAR", "1292 RAR\nBFs & Z"), outline=F, ylab="AIC",col=c("deepskyblue2", "brown"))
boxplot(RF~selZ, score, notch=T, names = c("282709\nnon RAR", "1292 RAR\nBFs & Z"), outline=F, ylab="normalized RF",col=c("deepskyblue2", "brown"))
dev.off()

#write.table(score,"score_gene_trees.txt", sep="\t", quote=F, row.names = F)

#ggplot version

p1 <- score %>% ggplot(aes(y=AIC, x=factor(sel), fill=factor(sel))) + geom_violin(scale="width", adjust=0.5, draw_quantiles=c(0.5), show.legend = FALSE) + theme_classic(base_size=14) + labs(x="Ratite Accelerated (Bayes Factor)?") + scale_x_discrete(labels=c("No", "Yes"))
p2 <- score %>% ggplot(aes(y=AIC, x=factor(selZ), fill=factor(selZ))) + geom_violin(scale="width", adjust=0.5, draw_quantiles=c(0.5), show.legend = FALSE) + theme_classic(base_size=14) + labs(x="Ratite Accelerated (Bayes Factor + Z)?") + scale_x_discrete(labels=c("No", "Yes"))
p3 <- score %>% ggplot(aes(y=RF, x=factor(sel), fill=factor(sel))) + geom_violin(scale="width", adjust=3, draw_quantiles=c(0.5), show.legend = FALSE) + theme_classic(base_size=14) + labs(x="Ratite Accelerated (Bayes Factor)?") + scale_x_discrete(labels=c("No", "Yes"))
p4 <- score %>% ggplot(aes(y=RF, x=factor(selZ), fill=factor(selZ))) + geom_violin(scale="width", adjust=3, draw_quantiles=c(0.5), show.legend = FALSE) + theme_classic(base_size=14) + labs(x="Ratite Accelerated (Bayes Factor + Z)?") + scale_x_discrete(labels=c("No", "Yes"))
                 
library(gridExtra)
s14<-grid.arrange(p1, p2, p3, p4, nrow=2)

ggsave("~/Projects/birds/ratite_compgen/manuscript/ScienceSubmissionRev1/FullDraftDec11/FigS14-Dec11.pdf", s14)
