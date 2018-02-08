#this should work from start to finish to regenerate the candidates that we refer to in the paper.

library(dplyr)
setwd("~/Desktop/candidates_summer2017/")

#load(file = ".RData")
#load atac
#atac <- read.table("~/Desktop/candidates_summer2017/tissue_specific_R_May31.txt", header=T, stringsAsFactors = F) %>% tbl_df
atac <- read.table("~/Desktop/tissue_specific_Nov9_update.txt", header=T, stringsAsFactors = F) %>% tbl_df

#get phylop and phyloacc

source("~/Desktop/candidates_summer2017/analyze_phyloP.R")
source("~/Desktop/candidates_summer2017/analyze_phyloAcc.R")

#merge and clean up
atac <- atac %>% select(CNEE:Strn9_strict, St16_whole_H3K27ac:TotalRepressed)
cnee <- full_join(phyloacc, atac, by=c("cnee" = "CNEE"))
cnee <- cnee %>% mutate(phylop.ratite = cnee %in% allRatite$name[allRatite$qval <= 0.05])
cnee <- cnee %>% mutate(phylop.kiwi = cnee %in% kiwi$name[kiwi$qval <= 0.05]) %>%
  mutate(phylop.moa = cnee %in% anoDid$name[anoDid$qval <= 0.05]) %>%
  mutate(phylop.cas = cnee %in% casuar$name[casuar$qval <= 0.05]) %>%
  mutate(phylop.ost = cnee %in% strCam$name[strCam$qval <= 0.05]) %>%
  mutate(phylop.rhea = cnee %in% rhea$name[rhea$qval <= 0.05]) %>%
  mutate(phylop.tin = cnee %in% tinamou$name[tinamou$qval <= 0.05]) %>%
  mutate(phylop.rcount = as.numeric(phylop.cas) + as.numeric(phylop.kiwi) + as.numeric(phylop.rhea) + as.numeric(phylop.moa) + as.numeric(phylop.ost))


#define properties - what we care about
#the .1 is if the bf is greater than 10
#the phylop.ratite is 1 if accelerated in ratites
#ratite_spec.1 need bf2 >1 and probability of non-ratite loss is low
#then we ask that it is not accelerated in tinamous.
#ratite_conv.1 requires ratite_loss_cons.prob is posterior estimate for number of independent losses. (needs to be greater than 2)
#ratite_conv.2 is stricter as it requires 2 losses in phylop as well

cnee <- cnee %>% mutate(ratite_accel.1 = bf1 >= 10 & ratite_loss.prob >= 1, ratite_accel.2 = ratite_accel.1 & phylop.ratite) %>%
  mutate(ratite_spec.1 = bf2 > 1 & nonratite_loss.prob < 1 & nonratite_loss.mat == 0, ratite_spec.2 = ratite_spec.1 & !phylop.tin) %>%
  mutate(ratite_conv.1 = ratite_loss_cons.prob > 2, ratite_conv.2 = ratite_conv.1 & ratite_loss_cons_min.mat >= 2 & phylop.rcount > 1)

#we need lots of new CNEE lists
#need start/stop for the cnees - 
cneeBED <- read.table("~/Dropbox/Doctorate/atac/MayJune_analysis/final_cnees_long.bed") %>% tbl_df %>% rename(cnee = V4,chromosome = V1,start=V2,stop=V3)
cnee <- full_join(cnee,cneeBED) %>% mutate(length=stop-start)

#ratite_accel.1 & ratite_spec.1 (bayesian model ratite accelerated, not exclusively convergent)
ratite.bayes.accel <- cnee %>% filter(ratite_accel.1 == T, ratite_spec.1 == T) %>% select(chromosome,start,stop)
write.table(x = ratite.bayes.accel,file = "oct27-ratite_bayes_accel.bed",col.names = F,row.names = F,quote = F,sep = "\t")

#ratite_accel.2 & ratite_spec.2 (subset of bayesian model also supported by phylop)
ratite.bayes.phylop.accel <- cnee %>% filter(ratite_accel.2 == T, ratite_spec.2 == T) %>% select(chromosome,start,stop)
write.table(x = ratite.bayes.phylop.accel,file = "oct27-ratite_bayes_phylo_accel.bed",col.names = F,row.names = F,quote = F,sep = "\t")

#ratite_accel.1 & ratite_spec.1 & ratite_conv.1 (subset of 1st list that are convergent)
ratite.bayes.conv <- cnee %>% filter(ratite_accel.1 == T, ratite_spec.1 == T,ratite_conv.1 == T) %>% select(chromosome,start,stop)
write.table(x = ratite.bayes.conv,file = "oct27-ratite_bayes_conv.bed",col.names = F,row.names = F,quote = F,sep = "\t")

#ratite_accel.2 & ratite_spec.2 & ratite_conv.2 (subset of 2nd list that are convergent)
ratite.bayes.phylop.conv <- cnee %>% filter(ratite_accel.2 == T, ratite_spec.2 == T,ratite_conv.2 == T) %>% select(chromosome,start,stop)
write.table(x = ratite.bayes.phylop.conv,file = "oct27-ratite_bayes_phylo_conv.bed",col.names = F,row.names = F,quote = F,sep = "\t")

#ratite_accel.1 & ratite_spec.1 & !ratite_conv.1 (subset of 1st list that are not convergent)
ratite.bayes.notconv <- cnee %>% filter(ratite_accel.1 == T, ratite_spec.1 == T, ratite_conv.1 == F) %>% select(chromosome,start,stop)
write.table(x = ratite.bayes.notconv,file = "oct27-ratite_bayes_notconv.bed",col.names = F,row.names = F,quote = F,sep = "\t")

#ratite_accel.2 & ratite_spec.2 & !ratite_conv.2 (subset of 2nd list that are not convergent)
ratite.bayes.phylop.notconv <- cnee %>% filter(ratite_accel.2 == T, ratite_spec.2 == T,ratite_conv.2 == F) %>% select(chromosome,start,stop)
write.table(x = ratite.bayes.phylop.notconv,file = "oct27-ratite_bayes_phylo_notconv.bed",col.names = F,row.names = F,quote = F,sep = "\t")

#check if CNEE length is available in this dataset.  pull 1000 at random (but NOT in any convergent or accel category above - 0 in that class).
#the minimum length should be the mean of the accelerated length
#run the enrichments on these as well.  
cnee <- cnee %>% mutate(bayes.accel=(ratite_accel.1 == T & ratite_spec.1 == T)) %>% mutate(bayes.accel2=(ratite_accel.2 == T & ratite_spec.2 == T))
table(cnee$bayes.accel,cnee$bayes.accel2)
test <- cnee %>% filter(bayes.accel==T) %>% select(length)
mean(test$length)
mean(cnee$length)
long.non.accel <- cnee %>% filter(bayes.accel==F,length>=195.5246) %>% select(chromosome,start,stop)
oneK.long.non.accel <- long.non.accel %>% sample_n(1000)

write.table(x = long.non.accel,file = "oct27-long_non_accel.bed",col.names = F,row.names = F,quote = F,sep = "\t")
write.table(x = oneK.long.non.accel,file = "oct27-long_1k_non_accel.bed",col.names = F,row.names = F,quote = F,sep = "\t")

cnee$length.code = cut(cnee$length, breaks=c(50,60,75,100,150,200,300,500,2265))
cnee %>% select(chromosome,start,stop,length.code) %>% write_tsv("cnee_cuts.bed", col_names = F)

non.accel <- cnee %>% filter(bayes.accel==F) %>% select(chromosome,start,stop,length) 
accel <- cnee %>% filter(bayes.accel==T) %>% select(chromosome,start,stop,length) 
accel

#load annotation
gene<-read.table("~/Desktop/candidates_summer2017/cnees.galgal4.annotation", header=F, stringsAsFactors = F) %>% tbl_df %>% rename(gene=V2)

cnee <- full_join(cnee, gene, by=c("cnee" = "V1"))

##FILTER FOR VENN DIAGRAM##

#all RARs used in the Venn Diagram: 760
final_cnee_cand_list_nofilter <- cnee %>% filter(ratite_accel.2, ratite_spec.2,ratite_conv.1 == T) %>% select(chromosome,start,stop)
#write.table(x = final_cnee_cand_list_nofilter,file = "nov14-ra2-rs2-rc1.bed",col.names = F,row.names = F,quote = F,sep = "\t")

#the 64 that fall into the middle of the Venn Diagram
final_cnee_cand_list <- cnee %>% filter(ratite_accel.2, ratite_spec.2,ratite_conv.1 == T,TotalActive >= 1, FL_strict > 0) %>% arrange(gene)
final_cnee_cand_list <- final_cnee_cand_list %>% mutate(FLspec = (TotalStrict == 1), broad = ((TotalStrict >= 5) & (TotalActive >= 5)))

# 26 that contain our 5 top candidates
final_cnee_cand_list %>% filter(FLspec | broad) %>% select(cnee, gene, TotalActive, TotalRepressed, ratite_loss_cons.mat, ratite_loss_cons_min.prob, bf1, bf2, phylop.ratite, phylop.rcount) %>% arrange(gene) %>% print.data.frame

