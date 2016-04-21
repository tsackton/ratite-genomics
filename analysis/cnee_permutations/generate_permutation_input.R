#generate permutation inputs
#currently requires analyze_cnees.R to be run first
#will eventually refactor to load data

perm.input<-accel.clade
perm.input[,ratiteQ:=min(c(Casuar, Kiwi, Ratite, Rhea, anoDid, strCam), na.rm=T), by=1:nrow(accel.clade)]
perm.input$ratite.acc = ifelse(perm.input$ratiteQ < 0.05 & perm.input$filtQ > 0.25, 1, 0)
perm.input.sub<-subset(perm.input, tinPres == 4 & nonRatiteLoss == 0, select=c(name, result.type, best_ens, ratite.acc, cas.acc, kiwi.acc, moa.acc, ostrich.acc, rhea.acc, cas.acc.2, kiwi.acc.2, ostrich.acc.2, rhea.acc.2))
perm.input.sub[,ratite.acc.2:=max(ratite.acc, cas.acc.2, kiwi.acc.2, ostrich.acc.2, rhea.acc.2), by=1:nrow(perm.input.sub)]

#clean up moa by setting to NA for type==wga
perm.input.sub[grepl("wga", result.type)]$moa.acc=NA

#the logic is that we are focusing on convergence in cases where we are (reasonably) confident that there is a ratite specific loss. We cannot rule out loss of constraint elsewhere, but by filtering for no CNEE losses outside of ratites in terms of sequence homology we limit the impact

#we will do one non-permutation based tests (compute FET for all pairwise sets of species for each analysis set), and a permutation based test where we keep track of the site pattern (01001) across ratites and compare to the empirical distribution. We will also do a by-gene permutation test, and a by-gene convergence test (the gene level tests use ratite.acc (1 or 2))

#all of this will be done on 12 input sets -- 6 result types x acc.1 | acc.2

#write out each type as a separate file
baseDir="/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/analysis/cnee_permutations/inputs/"

for (type in unique(perm.input.sub$result.type)) {
  write.table(perm.input.sub[result.type==type,c("best_ens", "ratite.acc", "cas.acc", "kiwi.acc", "moa.acc", "ostrich.acc", "rhea.acc"), with=F], file=paste0(baseDir, type, ".for_perm", ".1", sep=""), sep="\t", quote=F, row.names=F)
  write.table(perm.input.sub[result.type==type,c("best_ens", "ratite.acc.2", "cas.acc.2", "kiwi.acc.2", "ostrich.acc.2", "rhea.acc.2"), with=F], file=paste0(baseDir, type, ".for_perm", ".2", sep=""), sep="\t", quote=F, row.names=F)
}


