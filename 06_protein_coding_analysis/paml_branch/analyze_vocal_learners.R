source("setup_branch_runs.R")

##MAIN CODE TO RUN ANALYSIS AND WRITE OUT -- VOCAL LEARNERS

#command line options
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("Need to supply species to include and tree to analyze")
}

## VARIABLE DEFINITIONS ##

ts<-args[1]
tree_to_use<-args[2]
target_species <- c("melUnd", "calAnn", ts)
remove_species <- setdiff(c("pseHum", "ficAlb", "taeGut", "serCan", "geoFor", "corBra"),ts)
outfile = paste0(c("vl", ts, tree_to_use, "RER.out"), collapse="_")

compute_results(target_species, remove_species, tree_to_use) %>% write_tsv(path=outfile)


