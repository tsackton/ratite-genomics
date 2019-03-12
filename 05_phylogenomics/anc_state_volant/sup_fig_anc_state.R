#ancestral state reconstruction plotting - volant vs flightless tinamou ancestor

library(tidyverse)
#read in probs

vol<-read_csv("~/Projects/birds/ratite_compgen/ratite-genomics/05_phylogenomics/anc_state_volant/volant_state_probabilies_fine.csv")
names(vol) = c("loss", "gain", "palaeo", "eckrtm", "eckr", "eck", "mt")

sup_fig_ancstate <- vol %>% filter(gain <= 0.025) %>% ggplot(aes(y=gain, x=loss)) + geom_raster(aes(fill=mt), interpolate=TRUE) + scale_fill_gradient2(low="#0571b0", mid="#f7f7f7", high="#ca0020", limits=c(0,1), midpoint=0.5) + theme_gray(base_size = 14) + scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0,0))

#ggsave("~/Projects/birds/ratite_compgen/manuscript/ScienceSubmissionRev1/FullDraftDec11/ancState-Sup-Panel2.pdf", sup_fig_ancstate)

vol %>% filter(gain <= 0.025) %>% ggplot(aes(y=gain, x=loss)) + geom_raster(aes(fill=palaeo), interpolate=TRUE) + scale_fill_gradient2(low="#0571b0", mid="#f7f7f7", high="#ca0020", limits=c(0,1), midpoint=0.5) + theme_gray(base_size = 14) + scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0,0))

vol %>% filter(gain <= 0.025) %>% ggplot(aes(y=gain, x=loss)) + geom_raster(aes(fill=eckrtm), interpolate=TRUE) + scale_fill_gradient2(low="#0571b0", mid="#f7f7f7", high="#ca0020", limits=c(0,1), midpoint=0.5) + theme_gray(base_size = 14) + scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0,0))

vol %>% filter(gain <= 0.025) %>% ggplot(aes(y=gain, x=loss)) + geom_raster(aes(fill=eckr), interpolate=TRUE) + scale_fill_gradient2(low="#0571b0", mid="#f7f7f7", high="#ca0020", limits=c(0,1), midpoint=0.5) + theme_gray(base_size = 14) + scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0,0))

vol %>% filter(gain <= 0.025) %>% ggplot(aes(y=gain, x=loss)) + geom_raster(aes(fill=eck), interpolate=TRUE) + scale_fill_gradient2(low="#0571b0", mid="#f7f7f7", high="#ca0020", limits=c(0,1), midpoint=0.5) + theme_gray(base_size = 14) + scale_x_continuous(expand = c(0, 0)) + scale_y_reverse(expand = c(0,0))
