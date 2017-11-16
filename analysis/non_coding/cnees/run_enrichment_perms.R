#this code runs the permutations to test for GO, gene, and spatial enrichment in RARs and cRARs
library(tidyverse)
library(clusterProfiler)
library(org.Gg.eg.db)
library(parallel)

#set working directory
setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/")

#set cores and perms for whole run, can be individually modified later in code
CORES <- 5
PERMS <- 5

#load data
cnee <- read_tsv("cnees.tsv")

### GO ENRICHMENT HERE ###
#GO enrichment - four tests: accel .1s, accel .1s & conv .1s, accel .2s, accel. 1s & conv .1s & ratite_loss_cons_min.mat >= 2
#this code gets the real results
background <- cnee %>% filter(gene != ".") %>% dplyr::select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi)
set1 <- cnee %>% filter(gene != ".", ratite_accel.1, ratite_spec.1) %>% dplyr::select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi)
set2 <- cnee %>% filter(gene != ".", ratite_accel.2, ratite_spec.2) %>% dplyr::select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi)
set3 <- cnee %>% filter(gene != ".", ratite_accel.1, ratite_spec.1, ratite_conv.1) %>% dplyr::select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi)
set4 <- cnee %>% filter(gene != ".", ratite_accel.1, ratite_spec.1, ratite_conv.1, ratite_loss_cons_min.mat >= 2) %>% dplyr::select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi)
inputs <- list("set1" = set1, "set2" = set2, "set3" = set3, "set4" = set4)
calc_enrich <- function(targetset, background,ont) { enrichGO(targetset$ncbi,'org.Gg.eg.db',pvalueCutoff=1.5,qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, pAdjustMethod="none",universe=background,keytype="ENTREZID",ont=ont) }
bp_res_real <- lapply(inputs, calc_enrich, background=background$ncbi, ont="BP") %>% 
  lapply(slot, name="result") %>% 
  dplyr::bind_rows(.id = "set")
mf_res_real <-  lapply(inputs, calc_enrich, background=background$ncbi, ont="MF") %>% 
  lapply(slot, name="result") %>% 
  dplyr::bind_rows(.id = "set")
merged_mf_terms <- mf_res_real %>% dplyr::distinct(ID)
merged_bp_terms <- bp_res_real %>% dplyr::distinct(ID)

#get counts of CNEEs in each set
input_counts<-list("set1" = cnee %>% filter(gene != ".", ratite_accel.1, ratite_spec.1) %>% count %>% pull(n),
                   "set2" = cnee %>% filter(gene != ".", ratite_accel.2, ratite_spec.2) %>% count %>% pull(n),
                   "set3" = cnee %>% filter(gene != ".", ratite_accel.1, ratite_spec.1, ratite_conv.1) %>% count %>% pull(n),
                   "set4" = cnee %>% filter(gene != ".", ratite_accel.1, ratite_spec.1, ratite_conv.1, ratite_loss_cons_min.mat >= 2) %>% count  %>% pull(n))

#for each GO term in the merged go list, want to compute permutations: input is the merged term list and the counts

get_go_perm <- function(DF, samples, golist, ont) {
  rand <- DF %>% sample_n(samples) %>% filter(gene != ".") %>% 
    dplyr::select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
    distinct(ncbi)
  background <- DF %>% filter(gene != ".") %>% dplyr::select(gene) %>% separate(gene, into=c("ncbi", "sym"), sep=":") %>% distinct(ncbi)
  rand_go_bp <- calc_enrich(targetset=rand, background=background$ncbi, ont=ont)
  golist %>% left_join(rand_go_bp@result, by=c("ID" = "ID")) %>% separate(GeneRatio, into=c("target_in", "target_total")) %>% 
    separate(BgRatio, into=c("bg_in", "bg_total")) %>%
    mutate(newpval = ifelse(is.na(pvalue), 1, pvalue), logp.perm = -log10(newpval)) %>% 
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total), bg_frac = as.numeric(bg_in)/as.numeric(bg_total)) %>%
    dplyr::select(ID, logp.perm, target_frac, bg_frac) %>% arrange(ID)
}

#to do one permutation of the full set
get_one_perm_set <- function(perm, input, DF, golist, ont) {
  lapply(input, get_go_perm, DF=DF, golist=golist, ont=ont) %>%
  dplyr::bind_rows(.id="set")
}

#do permutations in parallel
para_cores <- CORES
num_perms <- PERMS

perm_bp <- mclapply(1:num_perms, get_one_perm_set, input=input_counts, DF=cnee, golist=merged_bp_terms, ont="BP", mc.cores=para_cores, mc.preschedule = TRUE) %>%
  dplyr::bind_rows(.id="perm")
perm_mf <- mclapply(1:num_perms, get_one_perm_set, input=input_counts, DF=cnee, golist=merged_mf_terms, ont="MF", mc.cores=para_cores, mc.preschedule = TRUE) %>%
  dplyr::bind_rows(.id="perm")

#write out
write_tsv(perm_bp, path="perm_bp_results.tsv")
write_tsv(perm_mf, path="perm_mf_results.tsv")

### GENE ENRICHMENT HERE ###
#for gene enrichment tests, the idea is to randomly permute each set and get count of CNEEs per gene



##JUNK BELOW HERE##
get_bp_perm(cnee, input_counts$set1, merged_bp_terms, "BP")

perm_bp_accel <- accel_go_bp@result %>% dplyr::select(ID) %>% arrange(ID)
bp_accel_count <- cnee %>% filter(gene != ".", ratite_accel.1, ratite_spec.1) %>% count
perm_bp_conv <- conv_go_bp@result %>% dplyr::select(ID) %>% arrange(ID)
bp_conv_count <- cnee %>% filter(gene != ".", ratite_accel.1, ratite_spec.1, ratite_conv.1) %>% count


for (n in 2:201) {
  perm_bp_accel[,n] = get_bp_perm(cnee, bp_accel_count$n, perm_bp_accel) %>% dplyr::select(logp.perm)
  perm_bp_conv[,n] = get_bp_perm(cnee, bp_conv_count$n, perm_bp_conv) %>% dplyr::select(logp.perm)
}

test <- left_join(accel_go_bp@result, perm_bp_logp, by=c("ID" = "GO"))
test %>% dplyr::select(ID, Description, pvalue, logp.perm:GO.9) %>% 
  group_by(ID) %>% filter(!is.na(GO.9)) %>% 
  mutate(logp.real = -log10(pvalue), max.perm = max(logp.perm:GO.9)) %>% arrange(desc(logp.real)) %>% filter(logp.real > max.perm+1) 

test2 <- left_join(conv_go_bp@result, perm_bp_logp_2, by=c("ID" = "GO"))
test2 %>% dplyr::select(ID, Description, pvalue, logp.perm:logp.perm.9) %>% 
  group_by(ID) %>% filter(!is.na(logp.perm.9)) %>% 
  mutate(logp.real = -log10(pvalue), max.perm = max(logp.perm:logp.perm.9)) %>% arrange(desc(logp.real)) %>% filter(logp.real > max.perm+1) 




#ready for analysis -- below code rought draft currently

##basic numbers
table(cnee$ratite_accel.1 & cnee$ratite_spec.1, cnee$ratite_accel.2 & cnee$ratite_spec.2)

cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count(num=floor(ratite_loss_cons_min.prob)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count(num=floor(ratite_loss_cons.prob)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count(num=floor(ratite_loss.prob)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count(num=floor(ratite_loss_cons_min.mat)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count(num=floor(ratite_loss_cons.mat)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count(num=floor(ratite_loss.mat)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))

cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% count(num=floor(ratite_loss_cons_min.prob)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% count(num=floor(ratite_loss_cons.prob)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% count(num=floor(ratite_loss.prob)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% count(num=floor(ratite_loss_cons_min.mat)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% count(num=floor(ratite_loss_cons.mat)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% count(num=floor(ratite_loss.mat)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))

#functional enrichment
#enrichment
library(clusterProfiler)
library(org.Gg.eg.db)


##PLOTS FOR TUFTS TALK##
library(ggthemes)
cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% ggplot(aes(x=floor(ratite_loss_cons.prob))) + geom_bar(fill="red") + coord_flip() + scale_x_reverse() + theme_classic() + theme(text=element_text(size=18))
#raw counts
cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count(floor(ratite_loss_cons_min.prob))


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

