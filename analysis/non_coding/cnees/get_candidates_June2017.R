library(tidyverse)

#load atac
atac <- read.table("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/atac/tissue_specific_R_May31.txt", header=T, stringsAsFactors = F) %>% tbl_df

#get phylop and phyloacc
source("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/analyze_phyloP.R")
source("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/analyze_phyloAcc.R")

#merge and clean up
atac <- atac %>% select(CNEE:Strn10_pool, St16_whole_H3K27ac:TotalRepressed)
#cnee <- full_join(phyloacc, atac, by=c("cnee" = "CNEE"))
cnee <- phyloacc %>% mutate(phylop.ratite = cnee %in% allRatite$name[allRatite$qval <= 0.05])
cnee <- cnee %>% mutate(phylop.kiwi = cnee %in% kiwi$name[kiwi$qval <= 0.05]) %>%
  mutate(phylop.moa = cnee %in% anoDid$name[anoDid$qval <= 0.05]) %>%
  mutate(phylop.cas = cnee %in% casuar$name[casuar$qval <= 0.05]) %>%
  mutate(phylop.ost = cnee %in% strCam$name[strCam$qval <= 0.05]) %>%
  mutate(phylop.rhea = cnee %in% rhea$name[rhea$qval <= 0.05]) %>%
  mutate(phylop.tin = cnee %in% tinamou$name[tinamou$qval <= 0.05]) %>%
  mutate(phylop.rcount = as.numeric(phylop.cas) + as.numeric(phylop.kiwi) + as.numeric(phylop.rhea) + as.numeric(phylop.moa) + as.numeric(phylop.ost))
  

#define properties

cnee <- cnee %>% mutate(ratite_accel.1 = bf1 >= 10 & ratite_loss.prob >= 1, ratite_accel.2 = ratite_accel.1 & phylop.ratite) %>%
  mutate(ratite_spec.1 = bf2 > 1 & nonratite_loss.prob < 1 & nonratite_loss.mat == 0, ratite_spec.2 = ratite_spec.1 & !phylop.tin) %>%
  mutate(ratite_conv.1 = ratite_loss_cons.prob > 2, ratite_conv.2 = ratite_conv.1 & ratite_loss_cons_min.mat >= 2 & phylop.rcount > 1)

#load annotation
gene<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/annotation/cnee/cnees.galgal4.annotation", header=F, stringsAsFactors = F) %>% tbl_df %>% rename(gene=V2)

cnee <- full_join(cnee, gene, by=c("cnee" = "V1"))

##FILTERS FOR CANDIDATES###

#FL specific putative enhancers -- broad
#cnee %>% filter(TotalStrict == 1, TotalActive >= 1, FL_strict > 0, ratite_accel.1, ratite_conv.1, ratite_spec.1) %>% select(cnee, FL_strict, TotalActive, TotalRepressed, ratite_loss_cons.mat, ratite_loss_cons.prob, bf1, bf2, gene, phylop.ratite, phylop.rcount) %>% arrange(ratite_loss_cons.mat, FL_strict) %>% print.data.frame

#FL specific putative enhancers -- narrow
#cnee %>% filter(TotalStrict == 1, TotalActive >= 1, FL_strict > 0, ratite_accel.2, ratite_conv.1, ratite_spec.2) %>% select(cnee, FL_strict, TotalActive, TotalRepressed, ratite_loss_cons.mat, ratite_loss_cons.prob, bf1, bf2, gene,phylop.ratite, phylop.rcount) %>% arrange(ratite_loss_cons.mat, FL_strict) %>% print.data.frame

#FL specific putative enhancers -- good compromise
#cnee %>% filter(TotalStrict == 1, TotalActive >= 1, FL_strict > 0, ratite_accel.2, ratite_spec.2, !phylop.tin, ratite_loss_cons_min.prob >= 1, ratite_loss_cons.mat >= 1) %>% select(cnee, gene, TotalActive, TotalRepressed, ratite_loss_cons.mat, ratite_loss_cons_min.prob, bf1, bf2, phylop.ratite, phylop.rcount) %>% arrange(gene) %>% print.data.frame

#Always on  putative enhancers -- good compromise
#cnee %>% filter(TotalStrict >= 5, TotalActive >= 5, ratite_accel.2, ratite_spec.2, !phylop.tin, ratite_loss_cons_min.prob > 1.5, ratite_loss_cons.mat > 1) %>% select(cnee, gene, TotalStrict, TotalActive, TotalRepressed, ratite_loss_cons.mat, ratite_loss_cons_min.prob, bf1, bf2, phylop.ratite, phylop.rcount) %>% arrange(gene) %>% print.data.frame


#drop chip-seq - 13 elements
#cnee %>% filter(AllBayes == 1, TotalStrict == 1, FL_strict > 0) %>% select(Chromosome:HL_strict, TotalActive, TotalRepressed) %>% arrange(Chromosome, Start) %>% print.data.frame

#relax convergent - 17 elements
#cnee %>% filter(set4_conv1 > 0 | set4_conv2 > 0 | set4_conv3 > 0 | set4_conv4 > 0, TotalStrict == 1, EverActive > 0, FL_strict > 0, TotalRepressed+TotalActive <= TotalActive | TotalActive > 0) %>% select(Chromosome:HL_strict, TotalActive, TotalRepressed) %>% arrange(Chromosome, Start) %>% print.data.frame

#relax FL-specificity - 43 elements
#cnee %>% filter(AllBayes == 1, TotalStrict > 0, EverActive > 0, FL_strict > 0, TotalRepressed+TotalActive <= TotalActive | TotalActive > 0) %>% select(Chromosome:HL_strict, TotalActive, TotalRepressed) %>% arrange(Chromosome, Start) %>% print.data.frame

#broadest set - 147 elements
#cnee %>% filter(set4_conv1 > 0 | set4_conv2 > 0 | set4_conv3 > 0 | set4_conv4 > 0, FL_strict > 0, EverActive > 0, TotalRepressed+TotalActive <= TotalActive | TotalActive > 0) %>% select(Chromosome:HL_strict, TotalActive, TotalRepressed) %>% arrange(Chromosome, Start)


##PLOTS FOR TUFTS TALK##
library(ggthemes)
cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% ggplot(aes(x=floor(ratite_loss_cons.prob))) + geom_bar(fill="red") + coord_flip() + scale_x_reverse() + theme_classic() + theme(text=element_text(size=18))
#raw counts
cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% count(floor(ratite_loss_cons.prob))

#neognath distribution of convergent losses
neo_losses <- postaccmat %>%
  mutate(taeGut.loss = (taeGut - taeGut.ficAlb)) %>%
  mutate(ficAlb.loss = (ficAlb - taeGut.pseHum)) %>%
  mutate(pseHum.loss = (pseHum - taeGut.pseHum)) %>%
  mutate(corBra.loss = (corBra - taeGut.corBra)) %>%
  mutate(melUnd.loss = (melUnd - taeGut.melUnd)) %>%
  mutate(falPer.loss = (falPer - taeGut.falPer)) %>%
  mutate(picPub.loss = (picPub - picPub.lepDis)) %>%
  mutate(lepDis.loss = (lepDis - picPub.lepDis)) %>%
  mutate(halLeu.loss = (halLeu - picPub.halLeu)) %>%
  mutate(aptFor.loss = (aptFor - aptFor.pygAde)) %>%
  mutate(pygAde.loss = (pygAde - aptFor.pygAde)) %>%
  mutate(fulGla.loss = (fulGla - aptFor.fulGla)) %>%
  mutate(nipNip.loss = (nipNip - aptFor.nipNip)) %>%
  mutate(balReg.loss = (balReg - balReg.chaVoc)) %>%
  mutate(chaVoc.loss = (chaVoc - balReg.chaVoc)) %>%
  mutate(calAnn.loss = (calAnn - calAnn.chaPel)) %>%
  mutate(chaPel.loss = (chaPel - calAnn.chaPel)) %>%
  mutate(cucCan.loss = (cucCan - calAnn.cucCan)) %>%
  mutate(colLiv.loss = (colLiv - colLiv.mesUni)) %>%
  mutate(mesUni.loss = (mesUni - colLiv.mesUni)) %>% select(taeGut.loss:mesUni.loss)

#sample 10000 tips in batches of 5 and get convergence counts
neo_rand_conv<-data.frame(rep=seq(1,10000), conv_ct_1=NA, conv_ct_2=NA, conv_ct_3=NA, conv_ct_4=NA, conv_ct_5=NA)
for (i in 1:10001) {
  targets <- sample(colnames(neo_losses), 5)
  conv <- neo_losses %>% select(targets) %>% mutate(conv = rowSums(.)) %>% select(conv)
  alllosses <- cbind(neo_losses %>% mutate(all_loss = rowSums(.)) %>% select(all_loss), conv)
  neo_rand_conv$conv_ct_1[i] <- sum(alllosses$conv == 1 & alllosses$all_loss == alllosses$conv) 
  neo_rand_conv$conv_ct_2[i] <- sum(alllosses$conv == 2 & alllosses$all_loss == alllosses$conv)
  neo_rand_conv$conv_ct_3[i] <- sum(alllosses$conv == 3 & alllosses$all_loss == alllosses$conv)
  neo_rand_conv$conv_ct_4[i] <- sum(alllosses$conv == 4 & alllosses$all_loss == alllosses$conv)
  neo_rand_conv$conv_ct_5[i] <- sum(alllosses$conv == 5 & alllosses$all_loss == alllosses$conv)
}

#real ratite data
sum(cnee$ratite_loss_cons.mat == 1 & cnee$tin_loss.mat == 0 & cnee$neo_loss.mat == 0)
sum(cnee$ratite_loss_cons.mat == 2 & cnee$tin_loss.mat == 0 & cnee$neo_loss.mat == 0)
sum(cnee$ratite_loss_cons.mat == 3 & cnee$tin_loss.mat == 0 & cnee$neo_loss.mat == 0)
sum(cnee$ratite_loss_cons.mat == 4 & cnee$tin_loss.mat == 0 & cnee$neo_loss.mat == 0)
sum(cnee$ratite_loss_cons.mat == 5 & cnee$tin_loss.mat == 0 & cnee$neo_loss.mat == 0)

#edit neo rand conv
neo_rand_conv <- neo_rand_conv %>% mutate(conv_tot = conv_ct_2 + conv_ct_3 + conv_ct_4 + conv_ct_5, total_accel = conv_tot + conv_ct_1)
neo_rand_conv <- neo_rand_conv %>% mutate(conv_prop = conv_tot / total_accel)

neo_rand_conv %>% ggplot(aes(x=conv_ct_2)) + geom_histogram(binwidth = 150)
neo_rand_conv %>% ggplot(aes(x=conv_ct_3)) + geom_histogram(binwidth = 10)
neo_rand_conv %>% ggplot(aes(x=conv_ct_4)) + geom_histogram(binwidth = 2)
neo_rand_conv %>% ggplot(aes(x=conv_ct_5)) + geom_histogram(binwidth = 1)

neo_rand_conv %>% ggplot(aes(x=conv_prop)) + geom_histogram(binwidth=0.01, fill="steelblue") + geom_vline(xintercept = (6+118+777+3578)/((6+118+777+3578)+18218), col="red") + theme_classic() + theme(text=element_text(size=18))
neo_rand_conv %>% 
  mutate(conv_prop2 = (conv_ct_3 + conv_ct_4 + conv_ct_5) / (conv_ct_1 + conv_ct_3 + conv_ct_4 + conv_ct_5)) %>% 
           ggplot(aes(x=conv_prop2)) + geom_histogram(binwidth=0.01, fill="steelblue") + geom_vline(xintercept = (6+118+777)/((6+118+777)+18218), col="red") + theme_classic() + theme(text=element_text(size=18))

cnee <- full_join(cnee, atac, by=c("cnee" = "CNEE"))

#testing enhancers generally
cnee <- cnee %>% mutate(ratite_set_1 = ratite_accel.1 & ratite_spec.1, ratite_set_2 = ratite_accel.1 & ratite_spec.1 & ratite_conv.1)
cnee %>% with(., table(TotalStrict > 0, ratite_set_1)) %>% fisher.test

cnee %>% with(., prop.table(table(TotalStrict > 0, ratite_set_1 & !ratite_set_2),2))
cnee %>% with(., prop.table(table(TotalStrict > 0, ratite_set_2),2))

cnee %>% with(., table(TotalStrict > 0 & HL_strict > 0, ratite_set_1)) %>% fisher.test
cnee %>% with(., table(TotalStrict > 0 & Keel10_strict > 0, ratite_set_1)) %>% fisher.test
cnee %>% with(., table(TotalStrict > 0 & Keel9_strict > 0, ratite_set_1)) %>% fisher.test
cnee %>% with(., table(TotalStrict > 0 & Strn10_strict > 0, ratite_set_1)) %>% fisher.test
cnee %>% with(., table(TotalStrict > 0 & Strn9_strict > 0, ratite_set_1)) %>% fisher.test
cnee %>% with(., table(TotalStrict > 0 & Pec10_strict > 0, ratite_set_1)) %>% fisher.test
cnee %>% with(., table(TotalStrict > 0 & Pec9_strict > 0, ratite_set_1)) %>% fisher.test

cnee %>% with(., table(TotalStrict == 1 & FL_strict > 0, ratite_set_1)) %>% fisher.test
cnee %>% with(., table(TotalStrict > 0 & HL_strict > 0, ratite_set_1)) %>% fisher.test

cnee %>% with(., table(TotalStrict == 1 & HL_strict > 0, ratite_set_1)) %>% fisher.test
cnee %>% with(., table(TotalStrict == 1 & Keel10_strict > 0, ratite_set_1)) %>% fisher.test
cnee %>% with(., table(TotalStrict == 1 & Keel9_strict > 0, ratite_set_1)) %>% fisher.test
cnee %>% with(., table(TotalStrict == 1 & Strn10_strict > 0, ratite_set_1)) %>% fisher.test
cnee %>% with(., table(TotalStrict == 1 & Strn9_strict > 0, ratite_set_1)) %>% fisher.test
cnee %>% with(., table(TotalStrict == 1 & Pec10_strict > 0, ratite_set_1)) %>% fisher.test
cnee %>% with(., table(TotalStrict == 1 & Pec9_strict > 0, ratite_set_1)) %>% fisher.test

#make gene lists for enrichment plots
cnee %>% filter(gene != ".") %>% select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi) %>% write_csv("all_cnee_genes.txt", col_names=FALSE)
cnee %>% filter(gene != ".", ratite_set_1) %>% select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi) %>% write_csv("ratite_accel_cnee_genes.txt", col_names=FALSE)
cnee %>% filter(gene != ".", ratite_set_2) %>% select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi) %>% write_csv("ratite_conv_cnee_genes.txt", col_names=FALSE)

cnee %>% filter(gene != ".") %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% filter(sym=="TBX5") %>% with(., table(ratite_set_2, FL_strict > 0))

#enrichment
biocLite("org.Gg.eg.db")
library(clusterProfiler)
library(org.Gg.eg.db)

ratite_accel_names <- cnee %>% filter(gene != ".", ratite_set_1) %>% select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi)
all_cnee_names <- cnee %>% filter(gene != ".") %>% select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi)
ratite_conv_names <- cnee %>% filter(gene != ".", ratite_set_2) %>% select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi) 

all_genes_k <- enrichKEGG(ratite_accel_names$ncbi,organism="gga",pvalueCutoff=0.05,pAdjustMethod="BH",qvalueCutoff=0.05,universe=all_cnee_names$ncbi,keyType="ncbi-geneid")
dotplot(all_genes_k)
enrichMap(all_genes_k)
plotGOgraph(all_genes_k)

accel_go_bp <- enrichGO(ratite_accel_names$ncbi,'org.Gg.eg.db',pvalueCutoff=0.05,pAdjustMethod="BH",qvalueCutoff=0.05,universe=all_cnee_names$ncbi,keytype="ENTREZID",ont="BP")
conv_go_bp <- enrichGO(ratite_conv_names$ncbi,'org.Gg.eg.db',pvalueCutoff=0.05,pAdjustMethod="BH",qvalueCutoff=0.05,universe=all_cnee_names$ncbi,keytype="ENTREZID",ont="BP")

accel_go_mf <- enrichGO(ratite_accel_names$ncbi,'org.Gg.eg.db',pvalueCutoff=0.05,pAdjustMethod="BH",qvalueCutoff=0.05,universe=all_cnee_names$ncbi,keytype="ENTREZID",ont="MF")
conv_go_mf <- enrichGO(ratite_conv_names$ncbi,'org.Gg.eg.db',pvalueCutoff=0.05,pAdjustMethod="BH",qvalueCutoff=0.05,universe=all_cnee_names$ncbi,keytype="ENTREZID",ont="MF")

accel_go_cc <- enrichGO(ratite_accel_names$ncbi,'org.Gg.eg.db',pvalueCutoff=0.05,pAdjustMethod="BH",qvalueCutoff=0.05,universe=all_cnee_names$ncbi,keytype="ENTREZID",ont="CC")
conv_go_cc <- enrichGO(ratite_conv_names$ncbi,'org.Gg.eg.db',pvalueCutoff=0.05,pAdjustMethod="BH",qvalueCutoff=0.05,universe=all_cnee_names$ncbi,keytype="ENTREZID",ont="CC")


dotplot(accel_go_bp)
dotplot(conv_go_bp)

dotplot(accel_go_mf)
dotplot(conv_go_mf)

