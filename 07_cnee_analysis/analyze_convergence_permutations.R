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
cnee %>% filter(ratite_spec.1, ratite_accel.1) %>% ggplot(aes(x=ratite_loss_cons.prob)) + 
  theme_classic() + geom_histogram(binwidth=0.1, center=2.05, closed="left") + labs(x="Posterior Number of Independent Ratite Accelerations") + 
  geom_vline(xintercept=2, col="red") + annotate("text", x=2.6, y=150, label = "839 convergent elements")
ggsave("~/Projects/birds/ratite_compgen/manuscript/ver5/ExtendedFigure7A.pdf")
cnee %>% filter(ratite_spec.1, ratite_accel.1) %>% ggplot(aes(x=ratite_loss_cons_min.prob)) + 
  theme_classic() + geom_histogram(binwidth=0.1, center=2.05, closed="left") + labs(x="Posterior Number of Independent Ratite Accelerations (Parsimony-restricted)") + 
  geom_vline(xintercept=2, col="red") + annotate("text", x=2.35, y=550, label = "521 convergent elements")
ggsave("~/Projects/birds/ratite_compgen/manuscript/ver5/ExtendedFigure7B.pdf")

cnee %>% filter(ratite_spec.1, ratite_accel.1) %>% ggplot(aes(x=ratite_loss_cons.mat)) + 
  theme_classic() + geom_histogram(binwidth=0.1, center=2.05, closed="left") + labs(x="Discrete Number of Independent Ratite Accelerations") + 
  geom_vline(xintercept=2, col="red") + annotate("text", x=2.65, y=1300, label = "586 convergent elements")
ggsave("~/Projects/birds/ratite_compgen/manuscript/ver5/ExtendedFigure7C.pdf")
cnee %>% filter(ratite_spec.1, ratite_accel.1) %>% ggplot(aes(x=ratite_loss_cons_min.mat)) + 
  theme_classic() + geom_histogram(binwidth=0.1, center=2.05, closed="left") + labs(x="Discrete Number of Independent Ratite Accelerations (Parsimony-restricted)") + 
  geom_vline(xintercept=2, col="red") + annotate("text", x=2.5, y=1400, label = "399 convergent elements")
ggsave("~/Projects/birds/ratite_compgen/manuscript/ver5/ExtendedFigure7D.pdf")

cnee %>% filter(ratite_spec.1, ratite_accel.1, ratite_loss_cons.prob >= 2) %>% count() 
cnee %>% filter(ratite_spec.1, ratite_accel.1, ratite_loss_cons.mat >= 2) %>% count() 
cnee %>% filter(ratite_spec.1, ratite_accel.1, ratite_loss_cons_min.prob >= 2) %>% count() 
cnee %>% filter(ratite_spec.1, ratite_accel.1, ratite_loss_cons_min.mat >= 2) %>% count() 


#null distribution Figure 3A
ggplot() + coord_cartesian(xlim=c(0,80)) + theme_classic() + 
  geom_density(adjust=3, fill="gray60",data=perm_conv, aes(x=count, ..count..)) + 
  geom_segment(aes(x=mean(obs_conv$conv_count), xend=mean(obs_conv$conv_count), y=50, yend=0), arrow=arrow(length=unit(0.5, "cm")), colour="red") +
  geom_jitter(data=obs_conv, aes(conv_count, y=1), colour="red", width=0, height=3, size=2) +
  annotate("text", x = mean(obs_conv$conv_count), y = 60, label = "Ratite mean")
ggsave("~/Projects/birds/ratite_compgen/manuscript/ver5/Figure3A.pdf")

