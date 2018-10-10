## DOING THIS UGLY BECAUSE DON'T HAVE THE TIME / ENERGY TO FIGURE OUT doParallel, ETC ##

source("setup_branch_runs.R")

##MAIN CODE TO RUN ANALYSIS AND WRITE OUT

#command line options
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("Need to supply species to include and tree to analyze")
}

## VARIABLE DEFINITIONS ##

ts<-args[1]
tree_to_use<-args[2]
target_species <- c("anoDid", "strCam", "nanAur", ts)
remove_species <- setdiff(c("rhePen", "rheAme", "aptHaa", "aptRow", "aptOwe", "droNov", "casCas"),ts)
outfile = paste0(c("ratite", ts, tree_to_use, "RER.out"), collapse="_")

compute_results(target_species, remove_species, tree_to_use) %>% write_tsv(path=outfile)


