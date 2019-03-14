#revised Oct 2018

library(tidyverse)
library(ape)

#functions
ctb <- function(x, cutoff = 0.90) {
  ifelse(x >= cutoff, 1, 0)
}

#load permutation results
wdir<-"~/Projects/birds/ratite_compgen/DRYAD/07_cnees/processed"

runs<-list()
for (whichset in c("extended", "original", "reduced")) {
    spec_patt<-glob2rx(paste0(whichset, "_perms_run*.tsv.gz"))
    files<-list.files(path=paste0(wdir, "/convperms"), pattern=spec_patt, full.names = TRUE)
    results<-list()
    for (file in files) {
      results[[file]] <- read_tsv(file) 
    }
    runs[[whichset]] <- bind_rows(results, .id="file")
}

perms<-bind_rows(runs, .id="set") %>% select(set, version, test, tips, count) %>% distinct(set, version, test, tips, .keep_all = TRUE)

#load real data
cnee_orig <- read_tsv(paste0(wdir, "/final_original_cnee.tsv.gz"))
cnee_red <- read_tsv(paste0(wdir, "/final_reduced_cnee.tsv.gz"))
cnee_ext <- read_tsv(paste0(wdir, "/final_extended_cnee.tsv.gz"))

phy_orig <- read.tree(paste0(wdir, "/original.phy"))
phy_red <- read.tree(paste0(wdir, "/reduced.phy"))
phy_ext <- read.tree(paste0(wdir, "/extended.phy"))

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
  
  return(tibble(count=count, tips=paste0(c(targets, cross_target), collapse="-")))
}

all_neo_orig <- cnee_orig %>% select(taeGut:anaPla) %>% names()
all_neo_ext <- cnee_ext %>% select(taeGut:nanBra,uriPel:anaPla) %>% names()
all_tin <- c("eudEle", "notPer", "tinGut", "cryCin")

#ratite dollo loop
ratite_gm_dollo <- list()
ratite_gl_dollo <- list()
for (third_sp in c("rheAme", "rhePen", "droNov", "casCas", "aptHaa", "aptOwe", "aptRow")) {
  test_targets <- c("strCam", "anoDid", third_sp)
  test_tips <- c(test_targets, all_neo_ext, all_tin)
  ratite_gm_dollo[[third_sp]] <- cnee_ext %>% filter(version=="gain") %>%
    conv_real(test_targets, test_tips)
  ratite_gl_dollo[[third_sp]] <- cnee_ext %>% filter(version=="gain_gap") %>%
    conv_real(test_targets, test_tips)
}

gm_dollo <- bind_rows(ratite_gm_dollo, .id="species")
gl_dollo <- bind_rows(ratite_gl_dollo, .id="species")


#new fig3A
perms %>% filter(version=="gain", set=="extended", test=="neo_conv_3") %>% ggplot(aes(count)) + 
  geom_histogram(binwidth = 1) + 
  coord_cartesian(xlim=c(0,50)) + theme_classic() + 
  geom_segment(aes(x=mean(gm_dollo$count), xend=mean(gm_dollo$count), y=50, yend=0), arrow=arrow(length=unit(0.5, "cm")), colour="red") +
  geom_jitter(data=gm_dollo, aes(count, y=4), colour="red", width=0, height=3, size=2) +
  ylab("Count") + xlab("Number of Convergently Accelerated Elements")


#alt versions
perms %>% filter(version=="gain_gap", set=="original", test=="neo_conv_3") %>% ggplot(aes(count)) + 
  geom_histogram(binwidth = 1) + 
  coord_cartesian(xlim=c(0,50)) + theme_classic() + 
  geom_segment(aes(x=mean(gl_dollo$count), xend=mean(gl_dollo$count), y=50, yend=0), arrow=arrow(length=unit(0.5, "cm")), colour="red") +
  geom_jitter(data=gl_dollo, aes(count, y=4), colour="red", width=0, height=3, size=2) +
  ylab("Count") + xlab("Number of Convergently Accelerated Elements")


#compute p-values
x=mean(gl_dollo$count)
perms %>% filter(version=="gain_gap", set=="original", test=="neo_conv_3") %>% 
  summarize(pval=(sum(count >= x)+1)/length(count))
x=median(gl_dollo$count)
perms %>% filter(version=="gain_gap", set=="original", test=="neo_conv_3") %>% 
  summarize(pval=(sum(count >= x)+1)/length(count))

x=mean(gm_dollo$count)
perms %>% filter(version=="gain", set=="original", test=="neo_conv_3") %>% 
  summarize(pval=(sum(count >= x)+1)/length(count))
x=median(gm_dollo$count)
perms %>% filter(version=="gain", set=="original", test=="neo_conv_3") %>% 
  summarize(pval=(sum(count >= x)+1)/length(count))

#numbers
mean(gl_dollo$count)
median(gl_dollo$count)
mean(gm_dollo$count)
median(gm_dollo$count)

#numbers from perms
perms %>% filter(test=="neo_conv_3") %>% with(., table(set, version))

perms %>% filter(test=="neo_conv_3", set=="original", version=="gain") %>% summarize(mean(count))

#counts / supplemental figures
#Fig S11 - plots of number of accelerated/convergent elements for gain, original dataset

cnee_gain_orig_conv <- cnee_orig %>% filter(dataset=="orig_v2_phyloAcc-gain", 
                     logBF1 >= 10, logBF2 > 1,
                     (it_pp_loss + ti_pp_loss) < 1, neo_tip_loss < 1) %>%
  mutate(floss_cl_bin = ctb(cd_pp_loss) + ctb(rh_pp_loss) + ctb(os_pp_loss) + ctb(ki_pp_loss) + ctb(mo_pp_loss),
         floss_cl_bin_dollo = ctb(ctb(cd_pp_loss) + ctb(rh_pp_loss) + ctb(ki_pp_loss)) + ctb(os_pp_loss) + ctb(mo_pp_loss)) %>%
  select(floss_cl_pp,floss_cl_pp_dollo,floss_cl_bin,floss_cl_bin_dollo,floss_sp_pp)

#different convergence metrics

#numbers

cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_bin > 1.8))
cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_bin_dollo > 1.8))
cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_pp > 1.8))
cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_pp_dollo > 1.8))

#percent
cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_bin > 1.8)/length(floss_cl_bin))
cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_bin_dollo > 1.8)/length(floss_cl_bin))
cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_pp > 1.8)/length(floss_cl_bin))
cnee_gain_orig_conv %>% summarize(count = sum(floss_cl_pp_dollo > 1.8)/length(floss_cl_bin))


#plots

cnee_gain_orig_conv %>% ggplot(aes(x=floss_cl_bin)) + geom_bar(fill="red") + theme_classic() + labs(x="Number of Accelerated Ratite Clades") + geom_vline(xintercept = 1.5)
ggsave("FigS11A.pdf")

cnee_gain_orig_conv %>% ggplot(aes(x=floss_cl_bin_dollo)) + geom_bar(fill="red") + theme_classic() + labs(x="Number of Accelerated Ratite Clades") + geom_vline(xintercept = 1.5)
ggsave("FigS11B.pdf")

cnee_gain_orig_conv %>% ggplot(aes(x=floss_cl_pp)) + geom_freqpoly(bins=50, col="red") + theme_classic() + labs(x="Posterior Estimate of Number of Independent Accelerations") + geom_vline(xintercept = 1.8)
ggsave("FigS11C.pdf")

cnee_gain_orig_conv %>% ggplot(aes(x=floss_cl_pp_dollo)) + geom_freqpoly(bins=50, col="red") + theme_classic() + labs(x="Posterior Estimate of Number of Independent Accelerations") + geom_vline(xintercept = 1.8)
ggsave("FigS11D.pdf")

#Fig S12 - consistency across models

cnee_orig %>% filter(logBF1 >= 10, logBF2 > 1, (it_pp_loss + ti_pp_loss) < 1) %>%
  ggplot(aes(floss_cl_pp, col=version)) + geom_freqpoly(bins=50, size=1.5) + theme_classic() +
  scale_color_brewer(palette = "Spectral")
ggsave("FigS12.pdf")


#Fig S13 - using gain version, compare original, extended, reduced

gain_comp_bf <- cnee_orig %>% filter(version=="gain") %>% select(cnee, logBF1, logBF2) %>% 
  inner_join(cnee_red %>% filter(version == "gain") %>% select(cnee, logBF1, logBF2), by=c("cnee" = "cnee"))

base_bf_plot <- gain_comp_bf %>% ggplot(aes(x=logBF1.x, y=logBF1.y)) + geom_hex(binwidth=c(3,3)) + theme_classic() + geom_hline(yintercept = 10, col="red") + geom_vline(xintercept = 10, col="red") + labs(x = "logBF1 (default dataset)", y="logBF1 (reduced dataset)") 

base_bf_plot + annotate("text", label = "Sig. in reduced, \n not in default: 404", x=-30, y=100) +
  annotate("text", label = "Sig. in both: 3379", x=200, y=100) +
  annotate("text", label = "Sig. in default, \n not in reduced: 1175", x=200, y=-25)
ggsave("FigS13A.pdf")

table(gain_comp_bf$logBF1.x >= 10, gain_comp_bf$logBF1.y >= 10)

gain_comp_pp <- cnee_orig %>% filter(version=="gain") %>% select(cnee, floss_cl_pp, logBF1, logBF2) %>% 
  full_join(cnee_red %>% filter(version == "gain") %>% select(cnee, floss_cl_pp, logBF1, logBF2), by=c("cnee" = "cnee")) %>% mutate(original = replace_na(floss_cl_pp.x, 0), reduced = replace_na(floss_cl_pp.y,0)) %>% 
  filter(logBF1.x >= 10 & logBF2.x >= 1 | logBF1.y >= 10 & logBF2.y >= 1)

gain_comp_pp %>% filter() %>% ggplot(aes(x=original, y=reduced)) + geom_point(alpha=0.2) + theme_classic()
ggsave("FigS13B.pdf")

