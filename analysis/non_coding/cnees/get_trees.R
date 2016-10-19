#load packages
library(ape)

#code to process CNEE branch length information
treeDir<-c("/Volumes/LaCie/Projects/Current/ratites/final/cnee_branchlen")
mods<-c("neutmod_ver1", "neutmod_ver2", "neutmod_ver3")
trees<-list()

for (neut in mods) {
  alltrees<-list.files(path=paste0(treeDir, "/", neut, collapse=""), recursive=T, full.names=T, pattern="RAxML_result.*")
  for (treefile in alltrees) {
    intree<-read.tree(treefile)
    treename<-sub(".*RAxML_result\\.", "", treefile, perl=T)
    treename<-sub("_ver\\d+", "", treename, perl=T)
    trees[[neut]][[treename]]<-intree
  }
}

save(trees, file=paste0(treeDir, "/processed_trees.Rlist"))
save(alltrees, file=paste0(treeDir, "/alltrees.Rlist"))

