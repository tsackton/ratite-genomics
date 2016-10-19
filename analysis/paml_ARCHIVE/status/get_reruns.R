
#setup
library(plyr)
setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/paml/status/")


hoglist<-read.table("all_hogs", header=F)
length(hoglist$V1)
names(hoglist)="hog"
models = c(0,1,2,7,8)
trees = c("True", "False")
modelset = data.frame(model_num = rep(models, length(trees)), species_tree = rep(trees, each=length(models)))

hog.status = data.frame(hog = rep(hoglist$hog, each=10), model_num=rep(modelset$model_num, length(hoglist)), species_tree=rep(modelset$species_tree, length(hoglist)))

site.runs<-read.table("site_parsed.out", header=T, sep="\t")

site.check<-merge(hog.status, site.runs[,c("hog", "species_tree", "model_num", "lnl")],all.x=T, all.y=F)
site.nummodels <- ddply(site.check, .(hog), summarize, num_models=length(unique(model_num[!is.na(lnl)])))

site.missing = subset(site.check, is.na(lnl))
site.missing = merge(site.missing, site.nummodels)
site.nummissing <- as.data.frame(table(site.missing$hog))
site.missing <- merge(site.missing, site.nummissing, by.x="hog", by.y="Var1")
#num models 5, freq 5 <- nothing missing
site.missing = subset(site.missing, !(num_models == 5 & Freq == 5))

#has species tree
#num models 5, Freq 1
#num models 4, Freq 2,3
#num models 3, Freq 4,5
#num models 2, Freq 6,7
#num models 1, Freq 8,9

#no species tree
#num models 4, Freq 6
#num models 3, Freq 7
#num models 2, Freq 8
#num models 1, Freq 9

site.missing<-subset(site.missing, !(species_tree == "True" & num_models == 4 & Freq == 6))
site.missing<-subset(site.missing, !(species_tree == "True" & num_models == 3 & Freq == 7))
site.missing<-subset(site.missing, !(species_tree == "True" & num_models == 2 & Freq == 8))

write.table(site.missing[,c("hog", "model_num", "species_tree")], file="site.reruns", sep="\t", quote=F, row.names=F, col.names=F)

#ancrec

ancrec.runs<-read.table("ancrec_parsed.out", header=T, sep="\t")
ancrec.check<-merge(hog.status, ancrec.runs[,c("hog", "species_tree", "treenum", "lnl")],all.x=T, all.y=F)
ancrec.check = subset(ancrec.check, model_num == 0)
nosp.hog = ancrec.check$hog[ancrec.check$treenum == 1 & ancrec.check$species_tree == "False"]
ancrec.missing = subset(ancrec.check, is.na(lnl) & !hog %in% nosp.hog)

write.table(ancrec.missing[,c("hog", "species_tree")], file="ancrec.reruns", sep="\t", quote=F, row.names=F, col.names=F)

#branch - site

#complicated to deal with trailing delimiter
branchsite.names<-c(as.character(read.table("branchsite_parsed.out", stringsAsFactors=F, nrows=1, sep="\t")),"drop")
branchsite.runs<-read.table("branchsite_parsed.out", skip=1, col.names=branchsite.names, sep="\t", comment.char="", header=F, stringsAsFactors=F, flush=T, check.names=F, colClasses=c(rep(NA,20), "NULL"))

branchsite.runs$foreground_species = gsub("([A-Za-z])_[A-Za-z0-9\\-_\\.]+", "\\1", branchsite.runs$foreground_species, perl=T)
branchsite.runs$foreground_species = apply(branchsite.runs[,"foreground_species", drop=F], 1, function(x) paste(unique(sort(unlist(strsplit(x, ":")))), sep="", collapse=":"))


#clean up duplicated trees
branchsite.runs$key = paste(branchsite.runs$model, branchsite.runs$hog, branchsite.runs$foreground_species, branchsite.runs$species_tree, sep=":")
branchsite.runs <- branchsite.runs[order(branchsite.runs$key, branchsite.runs$lnl, decreasing=T),]
branchsite.runs$duplicate = duplicated(branchsite.runs$key)
branchsite.clean = subset(branchsite.runs, duplicate==FALSE)
foreground_trees <- as.data.frame(table(branchsite.clean$foreground_species))

#extra classification of foreground trees
classify_string <- function(input) {
  invec <- unlist(strsplit(input, ":"))
  outvec <- character(0)
  kiwi <- character(0)
  cas <- character(0)
  ost <- character(0)
  rhea <- character(0)
  for (test in invec) {
    if (test == "aptRow" || test == "aptHaa" || test == "aptOwe") {
      kiwi <- "kiwi"  
    }
    if (test == "strCam") {
      ost <- "ostrich"      
    }
    if (test == "casCas" || test == "droNov") {
      cas <- "casuar"
    }
    if (test == "rheAme" || test == "rhePen") {
      rhea <- "rhea"
      
    }
  }
  outvec <- paste0(c(ost, kiwi, cas, rhea), collapse=":")
  return(outvec)  
}

foreground_trees$class = apply(foreground_trees[,c("Var1"), drop=F], 1, classify_string)
ddply(foreground_trees, .(class), summarize, sum(Freq))
foreground_trees$class[foreground_trees$class == "ostrich:kiwi:casuar:rhea" | foreground_trees$class == "ostrich:kiwi:casuar" | foreground_trees$class == "ostrich:kiwi:rhea" | foreground_trees$class == "ostrich:casuar:rhea" | foreground_trees$class == "ostrich:kiwi" | foreground_trees$class == "ostrich:rhea" | foreground_trees$class == "ostrich:casuar" | foreground_trees$class == "kiwi:casuar"] = "ratite"
foreground_trees$class[foreground_trees$class == "kiwi:casuar:rhea"] = NA


branchsite.clean <- merge(branchsite.clean, foreground_trees[,c("Var1", "class")], by.x="foreground_species", by.y="Var1", all.x=T, all.y=F)
branchsite.clean = subset(branchsite.clean, !is.na(branchsite.clean$class))

class = c("casuar", "kiwi", "ostrich", "ratite", "rhea")
trees = c("True", "False")
bs.modelset = data.frame(foreground_species = rep(class, length(trees)), species_tree = rep(trees, each=length(models)))
bs.status = data.frame(hog = rep(hoglist$hog, each=10), class=rep(bs.modelset$foreground_species, length(hoglist)), species_tree=rep(bs.modelset$species_tree, length(hoglist)))

bs.alt<-subset(branchsite.clean, model=="branchsite", select=c("hog", "class", "species_tree", "lnl"))
bs.null<-subset(branchsite.clean, model=="branchsitenull", select=c("hog", "class", "species_tree", "lnl"))

bs.alt<-merge(bs.alt, bs.status, all.x=F, all.y=T)
bs.null<-merge(bs.null, bs.status, all.x=F, all.y=T)

bs.alt.missing = subset(bs.alt, is.na(lnl) & !(species_tree == "True" & hog %in% nosp.hog))
bs.null.missing = subset(bs.null, is.na(lnl) & !(species_tree == "True" & hog %in% nosp.hog))

write.table(bs.alt.missing[,c("hog", "class", "species_tree")], file="bsalt.reruns", sep="\t", quote=F, row.names=F, col.names=F)
write.table(bs.null.missing[,c("hog", "class", "species_tree")], file="bsnull.reruns", sep="\t", quote=F, row.names=F, col.names=F)


#branch site initial analysis
bs.init<-subset(branchsite.clean, class=="ratite" & species_tree=="True" & !(hog %in% bs.alt.missing$hog[bs.alt.missing$class=="ratite" & bs.alt.missing$species_tree == "True"]) & !(hog %in% bs.null.missing$hog[bs.null.missing$class=="ratite" & bs.null.missing$species_tree == "True"])
bs.init<-bs.init[order(bs.init$hog, bs.init$model),c("hog", "class", "foreground_species", "model", "lnl")]

dup.hogs<-as.data.frame(table(bs.init$hog))[which(as.data.frame(table(bs.init$hog))[,2]!=2),][,1]
bs.init = subset(bs.init, !hog %in% dup.hogs)

library(reshape2)
bs.init$hog = as.character(bs.init$hog)
bs.wide<-dcast(bs.init, hog + class + foreground_species ~ model)
bs.wide$lrt = round(-2 * (bs.wide$branchsitenull - bs.wide$branchsite),4)
bs.wide$lrt[bs.wide$lrt<0]=0
bs.wide$pval = 1-pchisq(bs.wide$lrt, 1)
bs.wide$qval = p.adjust(bs.wide$pval, method="fdr")
