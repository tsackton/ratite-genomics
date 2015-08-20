#script to produce final version of guide tree for progressiveCactus with branch lengths
#based on UCE analysis of ratite genomes and the crocodile genome paper 4d tree

library(ape)
library(geiger)
library(phytools)

croc<-read.tree("crocpaper4d.tre") #4-fold degenerate tree from the crocodile genome paper
uce<-read.tree("ExaML_concatenated_best_tree_with_bootstraps_rerooted.tre") #preliminary UCE tree generated from our assemblies

##code to generate the UCE tree used here will be added in the future##

uce$tip.label=c("allMis", "anaPla", "melGal", "galGal", "mesUni", "colLiv", "chaPel", "calAnn", "cucCan", "chaVoc", "nipNip", "fulGla", "aptFor", "pygAde", "balReg", "picPub", "lepDis", "halLeu", "falPer", "melUnd", "corBra", "pseHum", "taeGut", "ficAlb", "casCas", "droNov", "aptHaa", "aptOwe", "aptRow", "rheAme", "rhePen", "notPer", "eudEle", "tinGut", "cryCin", "strCam")
croc$tip.label=sub("\\d+", "", croc$tip.label)

target<-read.tree("target.tre") #the topology derived from collapsing low confidence or conflicted nodes from Jarvis et al 2014

name.check(uce, data.names=target$tip.label)
name.check(croc, data.names=target$tip.label)

croc.clean<-drop.tip(croc, name.check(croc, data.names=target$tip.label)$tree_not_data)

uce.sub<-drop.tip(uce, name.check(uce, data.names=croc$tip.label)$tree_not_data)
uce.sub<-ladderize(uce.sub)
croc.sub<-drop.tip(croc, name.check(croc, data.names=uce$tip.label)$tree_not_data)
croc.sub<-ladderize(croc.sub)


plot(uce.sub, x.lim=0.6)
add.scale.bar()

plot(croc.sub, x.lim=0.6)
add.scale.bar()

br<-data.frame(croc=croc.sub$edge.length, uce=uce.sub$edge.length)
br<-br[br$croc<0.20,]
summary(lm(croc ~ 0 + uce, data=br))
plot(croc ~ uce, data=br, xlim=c(0,0.15), ylim=c(0,0.15))
abline(lm(croc ~ 0 + uce, data=br))
abline(a=0,b=1,lty="dashed",col="red")

#logic: use UCE tree for birds, croc tree for outgroups, but scale UCE branch length to be equivalent to croc branch lengths (i.e., multiple by 1.46583)

#first, collapse nodes to match target tree
uce.match<-uce
uce.match$node.label<-NULL
uce.match<-collapse.to.star(uce.match, fastMRCA(uce.match, "ficAlb", "pseHum"))


#
temp.edge<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "chaPel", "ficAlb"))]
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "chaPel", "ficAlb"))]=0
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "chaPel", "cucCan"))]<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "chaPel", "cucCan"))]+temp.edge
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "chaVoc", "ficAlb"))]<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "chaVoc", "ficAlb"))]+temp.edge

#
temp.edge<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "chaVoc", "ficAlb"))]
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "chaVoc", "ficAlb"))]=0
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "chaVoc", "pygAde"))]<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "chaVoc", "pygAde"))]+temp.edge
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "balReg", "ficAlb"))]<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "balReg", "ficAlb"))]+temp.edge

#
temp.edge<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "picPub", "halLeu"))]
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "picPub", "halLeu"))]=0
uce.match$edge.length[which(uce.match$edge[,2]==which(uce.match$tip.label=="halLeu"))]<-uce.match$edge.length[which(uce.match$edge[,2]==which(uce.match$tip.label=="halLeu"))]+temp.edge
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "picPub", "lepDis"))]<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "picPub", "lepDis"))]+temp.edge

#
temp.edge<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "balReg", "ficAlb"))]
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "balReg", "ficAlb"))]=0
uce.match$edge.length[which(uce.match$edge[,2]==which(uce.match$tip.label=="balReg"))]<-uce.match$edge.length[which(uce.match$edge[,2]==which(uce.match$tip.label=="balReg"))]+temp.edge
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "picPub", "ficAlb"))]<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "picPub", "ficAlb"))]+temp.edge

#
temp.edge<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "chaVoc", "pygAde"))]
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "chaVoc", "pygAde"))]=0
uce.match$edge.length[which(uce.match$edge[,2]==which(uce.match$tip.label=="chaVoc"))]<-uce.match$edge.length[which(uce.match$edge[,2]==which(uce.match$tip.label=="chaVoc"))]+temp.edge
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "pygAde", "nipNip"))]<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "pygAde", "nipNip"))]+temp.edge

#
temp.edge<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "chaPel", "cucCan"))]
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "chaPel", "cucCan"))]=0
uce.match$edge.length[which(uce.match$edge[,2]==which(uce.match$tip.label=="cucCan"))]<-uce.match$edge.length[which(uce.match$edge[,2]==which(uce.match$tip.label=="cucCan"))]+temp.edge
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "calAnn", "chaPel"))]<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "calAnn", "chaPel"))]+temp.edge

#
temp.edge<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "rhePen", "cryCin"))]
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "rhePen", "cryCin"))]=0
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "rhePen", "rheAme"))]<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "rhePen", "rheAme"))]+temp.edge
uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "notPer", "cryCin"))]<-uce.match$edge.length[which(uce.match$edge[,2]==fastMRCA(uce.match, "notPer", "cryCin"))]+temp.edge

plot(uce)
plot(uce.match)

uce.match<-di2multi(uce.match)
plot(uce.match)

#multiply branch lengths
uce.match$edge.length=uce.match$edge.length*1.46583

#write out trees
write.tree(croc.clean, file="croc.fixed.tre", digits=6)
write.tree(uce.match, file="uce.fixed.tre", digits=6)

#manually merge, using bird br from UCE tree and reptile branch lengths from croc tree; two crocs with no data use the branch lengths for allMis (for allsin) and croPor (for garhal)

#read back final manual tree
final<-read.tree(file = "final_tree.txt")
final.fix<-rotateConstr(final, target$tip.label)
write.tree(final.fix, file="final_tree_reorder.tre", digits=6)
