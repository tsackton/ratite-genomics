setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/")

#make list of all targets [tips and internal] to read in
targets<-read.table("phylop_target_list.txt", header=F, stringsAsFactors = F)
results<-list()
for (target in targets$V1) {
  infile=paste0(getwd(), "/phyloP_all_branches/", target, ".phyloP.out.gz")
  temp<-read.delim(infile)
  names(temp)[1]="cnee"
  temp$sp=target
  temp$qval=p.adjust(temp$pval, method="fdr")
  temp$padj=p.adjust(temp$pval, method="holm")
  results[[target]]=temp
}

phylop.fdr<-data.frame(cnee=results$anaPla$cnee)
for (target in targets$V1) {
  phylop.fdr[,target] = as.numeric(results[[target]]$qval<0.05)
}
row.names(phylop.fdr)<-phylop.fdr$cnee

phylop.holm<-data.frame(cnee=results$anaPla$cnee)
for (target in targets$V1) {
  phylop.holm[,target] = as.numeric(results[[target]]$padj<0.05)
}
row.names(phylop.holm)<-phylop.holm$cnee

#extract tips, repeat analysis
all.species<-read.delim("https://raw.githubusercontent.com/tsackton/ratite-genomics/master/alignment/species_list.tsv", stringsAsFactors = F)
vocal.tips<-c("taeGut", "ficAlb", "pseHum", "corBra", "melUnd", "calAnn")
ratite.tips<-c("aptRow", "aptOwe", "aptHaa", "casCas", "droNov", "rheAme", "rhePen", "anoDid", "strCam")
bird.tips<-c(all.species$Short.name[! all.species$Short.name %in% c("croPor", "gavGan", "anoCar", "chrPic","allMis","cheMyd", "allSin")], "anoDid")

tips.clean<-phylop.fdr[,bird.tips]

nonratite.losses = rowSums(tips.clean[,!colnames(tips.clean) %in% ratite.tips])
nonvl.losses = rowSums(tips.clean[,!colnames(tips.clean) %in% vocal.tips])
ratite.losses = rowSums(tips.clean[,ratite.tips])
all.losses<-data.frame(ratite=ratite.losses, vl=rowSums(tips.clean[,vocal.tips]), total=rowSums(tips.clean))

#load clade data


#permutation test


#three species test:
all_ratite_3sp=list()
for (rhea in c("rheAme", "rhePen")) {
  for (kiwi in c("aptHaa", "aptRow", "aptOwe")) {
    for (cas in c("droNov", "casCas")) {
      ratite_set = c(rhea, kiwi, cas, "strCam", "anoDid")
      ratite.3.set = combn(ratite_set, 3, simplif=F)
      all_ratite_3sp = c(all_ratite_3sp, ratite.3.set)
    }
  }
}
all_ratite_3sp = unique(all_ratite_3sp)

all_ratite_3sp_strict=list()
for (third_rat in c("rheAme", "rhePen","aptHaa", "aptRow", "aptOwe","droNov", "casCas")) {
  ratite_set = list(c(third_rat, "strCam", "anoDid"))
  all_ratite_3sp_strict = c(all_ratite_3sp_strict, ratite_set)
}


all_vl_3sp=list()
for (pass in c("taeGut", "ficAlb", "pseHum", "corBra")) {
  vl.3.set = list(c(pass, "melUnd", "calAnn"))
  all_vl_3sp=c(all_vl_3sp, vl.3.set)
}


#ratite counts
rat.3conv = numeric()
for (rat in all_ratite_3sp) {
  ratct<-data.frame(ratite=rowSums(tips.clean[,rat]),total=nonratite.losses)
  ratct=subset(ratct, total == 0 & ratite > 0)
  rat.3conv=c(rat.3conv,sum(ratct$ratite==3))
}

rat.3conv.strict=numeric()
for (rat in all_ratite_3sp_strict) {
  ratct<-data.frame(ratite=rowSums(tips.clean[,rat]),total=nonratite.losses)
  ratct=subset(ratct, total == 0 & ratite > 0)
  rat.3conv.strict=c(rat.3conv.strict,sum(ratct$ratite==3))
}

#vl counts
vl.3conv = numeric()
for (vl in all_vl_3sp) {
  ratct<-data.frame(ratite=rowSums(tips.clean[,vl]),total=nonvl.losses)
  ratct=subset(ratct, total == 0 & ratite > 0)
  vl.3conv=c(vl.3conv,sum(ratct$ratite==3))
}

#100 random sets
rand.3conv=numeric()
for (i in 0:100) {
  targets=sample(colnames(tips.clean))[1:3]
  ratct<-data.frame(ratite=rowSums(tips.clean[,targets]), total=rowSums(tips.clean))
  ratct=subset(ratct, ratite == total & ratite > 0)
  rand.3conv[i] = sum(ratct$ratite == 3)
}


