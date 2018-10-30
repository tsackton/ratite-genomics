#revised Oct 2018

library(tidyverse)
library(ape)

#set working dir
setwd("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/")

#load permutation results
wdir<-getwd()

runs<-list()
for (whichset in c("extended", "original", "reduced")) {
    spec_patt<-glob2rx(paste0(whichset, "_perms_run*.tsv"))
    files<-list.files(path=paste0(wdir, "/convperms"), pattern=spec_patt, full.names = TRUE)
    results<-list()
    for (file in files) {
      results[[file]] <- read_tsv(file) 
    }
    runs[[whichset]] <- bind_rows(results, .id="file")
}

perms<-bind_rows(runs, .id="set") %>% select(set, version, test, tips, count) %>% distinct(set, version, test, tips, .keep_all = TRUE)

#load real data
cnee_orig <- read_tsv("final_original_cnee.tsv.gz")
cnee_red <- read_tsv("final_reduced_cnee.tsv.gz")
cnee_ext <- read_tsv("final_extended_cnee.tsv.gz")

phy_orig <- read.tree("original.phy")
phy_red <- read.tree("reduced.phy")
phy_ext <- read.tree("extended.phy")

#compute observed convergences with same approach as permutations

conv_real <- function(DF, targets, tips, cutoff = 0.90) {
  count_acc_targets <- function(tip, prob_acc, cutoff, targets) {
    targets_selected <- prob_acc[tip %in% targets]
    return(sum(targets_selected > cutoff))
  }
  
  #compute 
  number<-length(targets)
  count<-DF %>% dplyr::select(cnee, one_of(tips)) %>% 
    gather(key = "tip", value = "prob_acc", -cnee) %>%
    group_by(cnee) %>%
    summarize(total_acc = sum(prob_acc > cutoff), target_acc = count_acc_targets(tip, prob_acc, cutoff, targets)) %>%
    dplyr::filter(target_acc == number & target_acc == total_acc) %>%
    tally() %>% pull(n)
  
  return(tibble(count=count, tips=paste0(targets, collapse="-")))
}

cross_real <- function(DF, targets, tips, cross_target, cutoff = 0.90) {
  count_acc_targets <- function(tip, prob_acc, cutoff, targets) {
    targets_selected <- prob_acc[tip %in% targets]
    return(sum(targets_selected > cutoff))
  }
  
  #compute 
  number<-length(targets)
  count<-DF %>% dplyr::select(cnee, one_of(tips)) %>% 
    gather(key = "tip", value = "prob_acc", -cnee) %>%
    group_by(cnee) %>%
    summarize(total_acc = sum(prob_acc > cutoff), 
              target_acc = count_acc_targets(tip, prob_acc, cutoff, targets), 
              cross_target_acc = count_acc_targets(tip, prob_acc, cutoff, cross_target)) %>%
    dplyr::filter(target_acc == number, target_acc == (total_acc-1), cross_target_acc == 1) %>%
    tally() %>% pull(n)
  
  return(tibble(count=count, tips=paste0(targets, cross_target, collapse="-")))
}

all_neo_orig <- cnee_orig %>% select(taeGut:anaPla) %>% names()
all_neo_ext <- cnee_ext %>% select(taeGut:nanBra,uriPel:anaPla) %>% names()
all_tin <- c("eudEle", "notPer", "tinGut", "cryCin")

#put this in a loop#

test_targets <- c("strCam")
test_tips <- c(test_targets, all_neo_ext, all_tin)

cnee_orig %>% filter(version=="gain") %>% conv_real(test_targets, test_tips)
cnee_ext %>% filter(version=="gain") %>% cross_real(test_targets, test_tips, "nanHar")

##### BELOW TO EDIT ###


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

