setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/paml_ancrec/")
library(data.table)

ancrec.muts<-fread("gunzip -c parsed_mutations.txt.gz")
names(ancrec.muts)<-c("hog", "tree", "branchpair", "ratite", "clean", "type", "mutclass")

ancrec.parsed<-fread("gunzip -c ancrec_parsed.out.gz")
ancrec.treekey<-ancrec.parsed[,c("hog", "treenum", "species_tree"), with=FALSE]
ancrec.muts<-merge(ancrec.muts, ancrec.treekey, by.x=c("hog", "tree"), by.y=c("hog", "treenum"))

#make gene tree subset 
ancrec.gt<-subset(ancrec.muts, clean=="clean" & species_tree==FALSE)

#make species tree subset
ancrec.st<-subset(ancrec.muts, clean=="clean" & species_tree==TRUE)

#clean=="clean" excludes branch pairs that are either sister to each other or where one is ancestral to the other#

#ANALYSIS FUNCTIONS, RUN WITH DIFFERENT INPUTS (gene tree, species tree) AND FILTERING OPTIONS
#SOME CODE REUSE BETWEEN FUNCTIONS COULD BE CLEANED UP BUT TOO LAZY##
#ALSO: NO ERROR CHECKING##

make_hist<-function(indata, freq = 10, subs = 100, file = "hist.pdf", type) {
  ancrec.clean <- indata
  #branchpair <-> ratite key
  branchpair.key<-unique(ancrec.clean[,c("branchpair", "ratite"), with=FALSE], by=c("branchpair", "ratite"))
  branchpair.ratite<-subset(branchpair.key, ratite == 2)
  branchpair.other<-subset(branchpair.key, ratite < 2)
  
  #branchpair <-> hog count
  branchpair.hog<-unique(ancrec.clean[,c("hog", "branchpair"), with=FALSE], by=c("hog", "branchpair"))
  branchpair.count<-as.data.frame(table(branchpair.hog$branchpair))
  branchpair.use<-as.character(branchpair.count$Var1[branchpair.count$Freq >= freq])
  
  if (type == "all") {
    ancrec.all.sum<-as.data.frame.matrix(table(ancrec.clean$branchpair, ancrec.clean$mutclass))
  } else if (type == "tips") {
    ancrec.clean.tips=subset(ancrec.clean, type=="tip-tip")
    ancrec.all.sum<-as.data.frame.matrix(table(ancrec.clean.tips$branchpair, ancrec.clean.tips$mutclass))
  }
  ancrec.all.sum$total = ancrec.all.sum$convergent + ancrec.all.sum$divergent + ancrec.all.sum$parallel + ancrec.all.sum$regular
  #filter
  ancrec.all.sum<-subset(ancrec.all.sum, row.names(ancrec.all.sum) %in% branchpair.use & total > subs)

  ancrec.all.sum$ratio = ancrec.all.sum$convergent / (ancrec.all.sum$divergent+1)
  ancrec.all.ratite = subset(ancrec.all.sum, rownames(ancrec.all.sum) %in% branchpair.ratite$branchpair)
  ancrec.all.other = subset(ancrec.all.sum, rownames(ancrec.all.sum) %in% branchpair.other$branchpair)  
  pdf(file=file)
  hist(ancrec.all.other$ratio, breaks=50, xlab="Convergent/Divergent Ratio", xlim=c(0,0.25), col="gray", main="")
  points(x=ancrec.all.ratite$ratio, y=rep(5, length(ancrec.all.ratite$ratio)), pch=16, col="red", cex=1)
  legend(x="right", legend=c("Ratite/Ratite branches"), pch=16, col="red", bty="n")
  dev.off()
}

make_hist(indata=ancrec.gt, freq=1, subs=25, file="aaconv_genetree_hist_freq1_subs25_allbranch.pdf", type="all")
make_hist(indata=ancrec.st, freq=1, subs=25, file="aaconv_sptree_hist_freq1_subs25_allbranch.pdf", type="all")
make_hist(indata=ancrec.gt, freq=1, subs=25, file="aaconv_genetree_hist_freq1_subs25_tipbranch.pdf", type="tips")
make_hist(indata=ancrec.st, freq=1, subs=25, file="aaconv_sptree_hist_freq1_subs25_tipbranch.pdf", type="tips")

make_hist(indata=ancrec.gt, freq=10, subs=100, file="aaconv_genetree_hist_freq10_subs100_allbranch.pdf", type="all")
make_hist(indata=ancrec.st, freq=10, subs=100, file="aaconv_sptree_hist_freq10_subs100_allbranch.pdf", type="all")
make_hist(indata=ancrec.gt, freq=10, subs=100, file="aaconv_genetree_hist_freq10_subs100_tipbranch.pdf", type="tips")
make_hist(indata=ancrec.st, freq=10, subs=100, file="aaconv_sptree_hist_freq10_subs100_tipbranch.pdf", type="tips")

make_hist(indata=ancrec.gt, freq=100, subs=100, file="aaconv_genetree_hist_freq100_subs100_allbranch.pdf", type="all")
make_hist(indata=ancrec.st, freq=100, subs=100, file="aaconv_sptree_hist_freq100_subs100_allbranch.pdf", type="all")
make_hist(indata=ancrec.gt, freq=100, subs=100, file="aaconv_genetree_hist_freq100_subs100_tipbranch.pdf", type="tips")
make_hist(indata=ancrec.st, freq=100, subs=100, file="aaconv_sptree_hist_freq100_subs100_tipbranch.pdf", type="tips")

get_conv_num<-function(indata) {
  ancrec.conv<-subset(indata, mutclass=="convergent" & type=="tip-tip")
  ancrec.conv.gene<-as.data.frame.matrix(table(ancrec.conv$hog, ancrec.conv$ratite==2))
  names(ancrec.conv.gene)=c("other", "ratite")
  return(table(ancrec.conv.gene$ratite))
}

get_conv_num(indata=ancrec.gt)
get_conv_num(indata=ancrec.st)

