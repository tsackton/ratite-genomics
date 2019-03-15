#Fig 13 (new) - ILS

library(tidyverse)

#read in discordance

ils<-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/04_ILS/gene_tree_discordance.txt") #sorry!

cnee<-read_tsv("~/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/final_original_cnee.tsv.gz") %>% filter(version == "gain")

ils<-full_join(ils, cnee, by=c("ID" = "cnee"))

ils <- ils %>% mutate(accelerated = ifelse(logBF1 >= 10 & logBF2 >= 1 & (ti_pp_loss + it_pp_loss) < 1, TRUE, FALSE))

p1 <- ils %>% ggplot(aes(y=AIC, x=factor(accelerated), fill=factor(accelerated))) + geom_violin(scale="width", adjust=0.5, draw_quantiles=c(0.5), show.legend = FALSE) + theme_classic(base_size=14) + labs(x="Ratite Accelerated?") + scale_x_discrete(labels=c("No", "Yes"))
p3 <- ils %>% ggplot(aes(y=RF, x=factor(accelerated), fill=factor(accelerated))) + geom_violin(scale="width", adjust=3, draw_quantiles=c(0.5), show.legend = FALSE) + theme_classic(base_size=14) + labs(x="Ratite Accelerated?") + scale_x_discrete(labels=c("No", "Yes"))

grid.arrange(p1,p3)
