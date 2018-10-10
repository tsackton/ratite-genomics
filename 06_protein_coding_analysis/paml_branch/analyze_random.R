## DOING THIS UGLY BECAUSE DON'T HAVE THE TIME / ENERGY TO FIGURE OUT doParallel, ETC ##

source("setup_branch_runs.R")

##MAIN CODE TO RUN ANALYSIS AND WRITE OUT

#command line options
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop("Need to supply number of tips to sample, tree to analyze, and iteration")
}

## VARIABLE DEFINITIONS ##

tree_to_use<-args[2]
tips_to_use <- setdiff(tree_to_use$tip.label,c("rhePen", "rheAme", "aptHaa", "aptRow", "aptOwe", "droNov", "casCas", "anoDid", "strCam"))
target_species <- get_random_targets(args[1], tips_to_use, trees[[tree_to_use]])
remove_species = NULL
outfile = paste0(c("random", args[1], args[3], tree_to_use, "RER.out"), collapse="_")

compute_results(target_species, remove_species, tree_to_use) %>% write_tsv(path=outfile)


