library(dplyr)
cnee <- read.table("~/Downloads/tissue_specific_R_May31.txt", header=T, stringsAsFactors = F) %>% tbl_df

##FILTERS FOR CANDIDATES###

#FL specific putative enhancers -- total of 6 elements to test
cnee %>% filter(AllBayes == 1, TotalStrict == 1, EverActive > 0, FL_strict > 0, TotalRepressed+TotalActive <= TotalActive | TotalActive > 0) %>% select(Chromosome:HL_strict, TotalActive, TotalRepressed) %>% arrange(FL_strict, HL_strict) %>% print.data.frame
fl_list <- cnee %>% filter(AllBayes == 1, TotalStrict == 1, EverActive > 0, FL_strict > 0, TotalRepressed+TotalActive <= TotalActive | TotalActive > 0) %>% select(CNEE)

#drop chip-seq - 13 elements
cnee %>% filter(AllBayes == 1, TotalStrict == 1, FL_strict > 0) %>% select(Chromosome:HL_strict, TotalActive, TotalRepressed) %>% arrange(Chromosome, Start) %>% print.data.frame

#relax convergent - 17 elements
cnee %>% filter(set4_conv1 > 0 | set4_conv2 > 0 | set4_conv3 > 0 | set4_conv4 > 0, TotalStrict == 1, EverActive > 0, FL_strict > 0, TotalRepressed+TotalActive <= TotalActive | TotalActive > 0) %>% select(Chromosome:HL_strict, TotalActive, TotalRepressed) %>% arrange(Chromosome, Start) %>% print.data.frame

#relax FL-specificity - 43 elements
cnee %>% filter(AllBayes == 1, TotalStrict > 0, EverActive > 0, FL_strict > 0, TotalRepressed+TotalActive <= TotalActive | TotalActive > 0) %>% select(Chromosome:HL_strict, TotalActive, TotalRepressed) %>% arrange(Chromosome, Start) %>% print.data.frame

#broadest set - 147 elements
cnee %>% filter(set4_conv1 > 0 | set4_conv2 > 0 | set4_conv3 > 0 | set4_conv4 > 0, FL_strict > 0, EverActive > 0, TotalRepressed+TotalActive <= TotalActive | TotalActive > 0) %>% select(Chromosome:HL_strict, TotalActive, TotalRepressed) %>% arrange(Chromosome, Start)
