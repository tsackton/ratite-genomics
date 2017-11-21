library("tidyverse")

#set working dir
setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/")

cnee <- read_tsv("cnees.tsv")
obs_conv <- read_tsv("~/Projects/birds/ratite_compgen/data/cnee_perms/obs_convergence.tsv")

perm_conv <-read_tsv("~/Projects/birds/ratite_compgen/data/cnee_perms/perm_convergence.tsv")

#remove row with all passerines as there are a lot of passerine-specific CNEEs losses so this is not convergence
perm_conv <- perm_conv %>% filter(count < 3000)

#add pval
obs_conv$epval = sapply(obs_conv$conv_count, function(x) (sum(perm_conv$count >= x)+1)/length(perm_conv$count))

#pvals for mean/median
x=mean(obs_conv$conv_count)
(sum(perm_conv$count >= x)+1)/length(perm_conv$count)
x=median(obs_conv$conv_count)
(sum(perm_conv$count >= x)+1)/length(perm_conv$count)

#null distribution stats
median(perm_conv$count)
mean(perm_conv$count)
length(perm_conv$count)

#counts - BF method
conv_counts<-list("set1" = cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% count %>% pull(n),
                   "set2" = cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count %>% pull(n),
                   "set3" = cnee %>% filter(ratite_accel.1, ratite_spec.1, ratite_conv.1) %>% count %>% pull(n),
                   "set4" = cnee %>% filter(ratite_accel.1, ratite_spec.1, ratite_conv.1, ratite_loss_cons_min.mat >= 2) %>% count  %>% pull(n))

#plots [extended figure XX]
cnee %>% filter(ratite_spec.1, ratite_accel.1) %>% ggplot(aes(x=ratite_loss_cons.prob)) + geom_density()
cnee %>% filter(ratite_spec.1, ratite_accel.1) %>% ggplot(aes(x=ratite_loss_cons_min.mat)) + geom_histogram()

cnee %>% filter(ratite_spec.1, ratite_accel.1, ratite_loss_cons_min.prob >= 2) %>% count() 
cnee %>% filter(ratite_spec.1, ratite_accel.1, ratite_loss_cons.mat >= 2) %>% count() 
cnee %>% filter(ratite_spec.1, ratite_accel.1, ratite_loss_cons_min.mat >= 2) %>% count() 


#null distribution Figure 3A
ggplot() + coord_cartesian(xlim=c(0,80)) + theme_classic() + 
  geom_density(adjust=3, fill="gray60",data=perm_conv, aes(x=count, ..count..)) + 
  geom_segment(aes(x=mean(obs_conv$conv_count), xend=mean(obs_conv$conv_count), y=50, yend=0), arrow=arrow(length=unit(0.5, "cm")), colour="red") +
  geom_jitter(data=obs_conv, aes(conv_count, y=1), colour="red", width=0, height=3, size=2) +
  annotate("text", x = mean(obs_conv$conv_count), y = 60, label = "Ratite mean")


