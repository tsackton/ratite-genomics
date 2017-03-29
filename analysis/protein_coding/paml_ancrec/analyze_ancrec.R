setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/paml_ancrec/")
library(data.table)
library(ape)

#read in files
ancrec.parsed<-fread("gunzip -c paml_M0_parsed.txt.gz")
ancrec.treekey<-ancrec.parsed[,c("hog", "treenum", "species_tree"), with=FALSE]
hog_info<-read.table("../all_hog_info.tsv", sep="\t", header=T)

hog_info$has_species_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree]
hog_info$has_gene_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree == F]

ancrec.muts<-fread("gunzip -c paml_M0_parsedmuts.txt.gz")
names(ancrec.muts)<-c("hog", "tree", "branchnames", "branchpair", "ratite", "clean", "type", "mutclass")

ancrec.muts.merge<-merge(ancrec.muts, ancrec.treekey, by.x=c("hog", "tree"), by.y=c("hog", "treenum"), all=T)

#get tree
fulltree<-read.tree("../final_neut_tree.nwk")
tipdist<-cophenetic.phylo(fulltree)

#recaculate missing after removing species not in tree
tips_to_keep = fulltree$tip.label[fulltree$tip.label %in% names(hog_info)]
hog_info$missing_ct_clean = apply(hog_info[,tips_to_keep], 1, function(x) sum(is.na(x)))

#merge
ancrec.muts.merge<-merge(ancrec.muts.merge, hog_info, by="hog")

#okay, let's start with a conservative set: species tree, no gene duplications, only tips
#redo filtering to remove species not included in tree prior to computing missing and dup_ct measures

ancrec.1<-subset(ancrec.muts.merge, species_tree == TRUE & dup_ct == 0 & clean=="clean" & type=="tip-tip" & missing_ct_clean<3, select=c("hog", "branchnames", "branchpair", "ratite", "mutclass"))
#get number of hogs
length(unique(ancrec.1$hog))

ancrec.1.sum<-as.data.frame.matrix(table(ancrec.1$branchpair, ancrec.1$mutclass))
ancrec.1.key<-unique(ancrec.1[,c("branchpair", "branchnames", "ratite"), with=FALSE], by=c("branchpair", "branchnames", "ratite"))

#add between clade and within clade info to key: neo-neo, paleo-paleo, neo-paleo
palaeo<-c("rheAme", "rhePen", "droNov", "casCas", "aptRow", "aptOwe", "aptHaa", "strCam", "tinGut", "cryCin", "notPer", "eudEle")
ancrec.1.key$branch1 = sub("-\\w+", "", ancrec.1.key$branchnames)
ancrec.1.key$branch2 = sub("\\w+-", "", ancrec.1.key$branchnames)
ancrec.1.key$dist = apply(ancrec.1.key, 1, function(x) if(x[4] %in% rownames(tipdist) & x[5] %in% rownames(tipdist)) { tipdist[x[4], x[5]] } else { return(NA) })

ancrec.1.key$clade = "neo-palaeo"
ancrec.1.key$clade[ancrec.1.key$branch1 %in% palaeo & ancrec.1.key$branch2 %in% palaeo] = "palaeo-palaeo"
ancrec.1.key$clade[!(ancrec.1.key$branch1 %in% palaeo) & !(ancrec.1.key$branch2 %in% palaeo)] = "neo-neo"

ancrec.1.sum<-merge(ancrec.1.sum, ancrec.1.key, by.x="row.names", by.y="branchpair")
ancrec.1.sum$color="black"
ancrec.1.sum$color[ancrec.1.sum$ratite==2]="red"
ancrec.1.sum$broadconv = ancrec.1.sum$parallel+ancrec.1.sum$convergent
ancrec.1.sum$broaddiv = ancrec.1.sum$divergent + ancrec.1.sum$regular
ancrec.1.sum$broadratio = ancrec.1.sum$broadconv/(ancrec.1.sum$broaddiv+1)
ancrec.1.sum$totalsubs = ancrec.1.sum$broadconv + ancrec.1.sum$broaddiv

#done prepping data for species tree analysis

#linear model
mod.full<-lm(broadconv ~ broaddiv*dist + as.factor(ratite==2), data=ancrec.1.sum)
mod.noratite<-lm(broadconv ~ broaddiv*dist, data=ancrec.1.sum)
anova(mod.noratite,mod.full)

#figure 2a

plot(broadconv ~ broaddiv, col=color, pch=16, data=ancrec.1.sum[ancrec.1.sum$dist < 0.15,], log="xy", xlab="# of Divergent Subs", ylab="# of Convergent Subs", bty="l")
nonratite.mod<-lm(broadconv ~ 0 + broaddiv, data=ancrec.1.sum[ancrec.1.sum$ratite < 2 & ancrec.1.sum$dist<0.15,])

plotpoints=data.frame(x=seq(1,10000,10), y=seq(1,10000,10)*coef(nonratite.mod))
points(plotpoints, type='l')

legend(x="topleft", legend=c("Ratite-ratite branch pair", "Other branch pairs"), pch=16, col=c("red", "black"), bty="n")

#specific genes
get_conv_num<-function(indata) {
  ancrec.conv<-subset(indata, (mutclass=="convergent" | mutclass=="parallel"))
  ancrec.all<-indata
  ancrec.all.gene<-as.data.frame.matrix(table(ancrec.all$hog, ancrec.all$ratite==2))
  names(ancrec.all.gene)=c("non", "rat")
  ancrec.all.gene$prop.ratite = ancrec.all.gene$rat/(ancrec.all.gene$non+ancrec.all.gene$rat)
  ancrec.conv.gene<-as.data.frame.matrix(table(ancrec.conv$hog, ancrec.conv$ratite==2))
  names(ancrec.conv.gene)=c("non", "rat")
  ancrec.test<-merge(ancrec.all.gene, ancrec.conv.gene, by="row.names")
  #conv is y
  ancrec.test2<-subset(ancrec.test, prop.ratite > 0)
  ancrec.test2$pval<-apply(ancrec.test2, 1, function(x) prop.test(as.numeric(x[6]),as.numeric(x[5])+as.numeric(x[6]),as.numeric(x[4]),alt="g")$p.value)
  return(ancrec.test2)
}

ratite.conv.genes<-get_conv_num(indata=ancrec.1)

##robustness

##STILL TO FINALIZE
#REDO?

