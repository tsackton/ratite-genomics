#analyze cons_scale results

setwd("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/")
csrows<-scan("cons_scale_ver2_rows", what="character")
cscols<-scan("cons_scale_ver2_header", what="character")
cons<-read.table("cons_scale_ver2_elem_Z.txt", stringsAsFactors=F, header=F, row.names=csrows, col.names=cscols)
head(cons)
ratite<-c("aptRow", "aptOwe", "aptHaa", "casCas", "droNov", "rheAme", "rhePen", "anoDid", "strCam", "aptHaa.aptOwe", "aptHaa.aptRow", "casCas.droNov", "rheAme.rhePen")
nonbird<-c("allMis", "allSin", "croPor", "gavGan", "chrPic", "cheMyd", "anoCar", "taeGut.aptHaa", "allMis.allSin", "croPor.gavGan", "allMis.croPor", "taeGut.allMis", "chrPic.cheMyd", "taeGut.chrPic", "taeGut.anoCar")
paleo<-c("cryCin", "tinGut", "eudEle", "notPer", "aptHaa.casCas", "aptHaa.rheAme", "cryCin.tinGut", "eudEle.notPer", "cryCin.eudEle", "cryCin.anoDid", "aptHaa.cryCin", "aptHaa.strCam")
penguin<-c("pygAde", "aptFor", "aptFor.pygAde")
vocal.tips<-c("taeGut", "ficAlb", "pseHum", "corBra", "melUnd", "calAnn")
ratite.tips<-c("aptRow", "aptOwe", "aptHaa", "casCas", "droNov", "rheAme", "rhePen", "anoDid", "strCam")
neognath<-colnames(cons)[!(colnames(cons) %in% nonbird | colnames(cons) %in% ratite | colnames(cons) %in% paleo | colnames(cons) %in% penguin)]
allbirds<-c(neognath,penguin,ratite,paleo)
cons$ratite<-apply(cons[,ratite], 1, function(x) sum(x != 1))
cons$nonbird<-apply(cons[,nonbird], 1, function(x) sum(x != 1))
cons$penguin<-apply(cons[,penguin], 1, function(x) sum(x != 1))
cons$paleo<-apply(cons[,paleo], 1, function(x) sum(x != 1))
cons$neognath<-apply(cons[,neognath], 1, function(x) sum(x != 1))

#for simplicity, let's analyze just tips
tips<-cons[,c(1:36)]

#proof of principle
#take a random sample of neognaths
random_set=data.frame(ct0=numeric(), ct1=numeric(), ct2=numeric(), ct3=numeric(), ct4=numeric(), ct5=numeric())
ratite_set=data.frame(ct0=numeric(), ct1=numeric(), ct2=numeric(), ct3=numeric(), ct4=numeric(), ct5=numeric())

for (i in 1:100) {
  target<-sample(neognath[!grepl(".", neognath, fixed=T)])[1:5]
  res<-data.frame(nontarget=apply(tips[,!(colnames(tips) %in% target)], 1, function(x) sum(x != 1)), target=factor(apply(tips[,target], 1, function(x) sum (x != 1)), levels=c(0,1,2,3,4,5)))
  random_set[i,]=table(res$target[res$nontarget==0])[1:6]
  target<-sample(ratite[!grepl(".", ratite, fixed=T)])[1:5]
  res<-data.frame(nontarget=apply(tips[,!(colnames(tips) %in% target)], 1, function(x) sum(x != 1)), target=factor(apply(tips[,target], 1, function(x) sum (x != 1)), levels=c(0,1,2,3,4,5)))
  ratite_set[i,]=table(res$target[res$nontarget==0])[1:6]
}

tips$total.loss = apply(tips, 1, function(x) sum(x != 1))
tips$ratite.loss = apply(tips[,ratite[!grepl(".", ratite, fixed=T)]], 1, function(x) sum(x != 1))
tips$vl.loss = apply(tips[,c("taeGut", "ficAlb", "pseHum", "corBra", "melUnd", "calAnn")], 1, function(x) sum (x != 1))

tips[tips$total.loss==9 & tips$ratite.loss == 9,]

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

tips.clean=tips[,c(1:36)]
tips.clean[tips.clean==1]=5
tips.clean[tips.clean<5]=1
tips.clean[tips.clean==5]=0

nonratite.losses = rowSums(tips.clean[,!colnames(tips.clean) %in% ratite.tips])
nonvl.losses = rowSums(tips.clean[,!colnames(tips.clean) %in% vocal.tips])
ratite.losses = rowSums(tips.clean[,ratite.tips])
all.losses<-data.frame(ratite=ratite.losses, vl=rowSums(tips.clean[,vocal.tips]), total=rowSums(tips.clean))
nonbird.losses<-rowSums((tips.clean[,nonbird]))

#candidates
ratite.cands<-rownames(all.losses[all.losses$ratite == all.losses$total & all.losses$ratite > 4,])
ratite.cands.loose<-rownames(all.losses[all.losses$ratite == all.losses$total & all.losses$ratite > 1,])
ratite.cands.strict<-rownames(all.losses[all.losses$ratite == all.losses$total & all.losses$ratite >= 8,])
ratite.cands.all<-rownames(all.losses[all.losses$ratite == all.losses$total & all.losses$ratite >= 1,])

#add info
ce.info<-read.csv("~/Dropbox/Work/EggCNEEs/October 2016/gg.final.annotated.csv")
names(ce.info)[2]="cnee"
ce.info$rar = ce.info$cnee %in% ratite.cands
ce.info$rar.s = ce.info$cnee %in% ratite.cands.strict
ce.info$rar.l = ce.info$cnee %in% ratite.cands.loose
ce.info$rar.a = ce.info$cnee %in% ratite.cands.all

write.table(unique(ce.info$best_ens[ce.info$rar.a == T]), file="rar.all.ensID", quote=F, sep="\t", row.names=F, col.names=F)
write.table(unique(ce.info$best_ens[ce.info$rar == T]), file="rar.med.ensID", quote=F, sep="\t", row.names=F, col.names=F)
write.table(unique(ce.info$best_ens[ce.info$rar.s == T]), file="rar.strict.ensID", quote=F, sep="\t", row.names=F, col.names=F)
write.table(unique(ce.info$best_ens[ce.info$rar.l == T]), file="rar.loose.ensID", quote=F, sep="\t", row.names=F, col.names=F)
write.table(unique(ce.info$best_ens), file="rar_background.ensID", quote=F, sep="\t", row.names=F, col.names=F)
table(ce.info$rar.a)

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

#10000 random sets
rand.3conv=numeric()
for (i in 0:10000) {
  targets=sample(colnames(tips.clean))[1:3]
  ratct<-data.frame(ratite=rowSums(tips.clean[,targets]), total=rowSums(tips.clean))
  ratct=subset(ratct, ratite == total & ratite > 0)
  rand.3conv[i] = sum(ratct$ratite == 3)
}


#save permutation
write.table(rand.3conv, file="rand.3conv.perms", quote=F, sep="\t", row.names=F, col.names=F)

#read permutation
rand.3conv<-read.table("rand.3conv.perms", sep="\t", stringsAsFactors = F, header=F)

#compare
plot(y=seq(10000,1,-1)/10000,x=sort(rand.3conv$V1), type="l", col="black", lwd=3, lty="dashed", xlim=c(0,100), las=1, xlab="Number Convergent", ylab="Proportion")
points(y=0.2, x=mean(rat.3conv), col="blue", pch=16, cex=1.5)
points(y=0.2, x=mean(rat.3conv.strict), col="blue", pch=17, cex=1.5)
points(y=0.2, x=mean(vl.3conv), col="firebrick", pch=16, cex=1.5)
legend("topright", legend=c("Vocal Learners", "Ratites", "Ratites [strict]"), pch=c(16,16,17), col=c("firebrick", "blue", "blue"))
text(x=mean(vl.3conv)-1, y=0.25, labels=paste("P=", sum(rand.3conv >= mean(vl.3conv))/10000), adj=c(0,1))
text(x=mean(rat.3conv.strict)-1, y=0.25, labels=paste("P=", sum(rand.3conv >= mean(rat.3conv.strict))/10000), adj=c(0.3,1))
text(x=mean(rat.3conv)-1, y=0.25, labels=paste("P=", sum(rand.3conv >= mean(rat.3conv))/10000), adj=c(0,1))

table(all.losses$ratite[all.losses$total == all.losses$ratite & all.losses$total > 0])
ratite.conv.all.cnees<-rownames(all.losses[all.losses$total == all.losses$ratite & all.losses$total > 0,])
tips.ratite.losses<-tips.clean[ratite.conv.all.cnees,c("anoDid", "strCam", "rheAme", "rhePen", "casCas", "droNov", "aptHaa", "aptRow", "aptOwe")]
tips.ratite.losses$rhea = pmax(tips.ratite.losses$rheAme,tips.ratite.losses$rhePen)
tips.ratite.losses$kiwi = pmax(tips.ratite.losses$aptOwe,tips.ratite.losses$aptRow,tips.ratite.losses$aptHaa)
tips.ratite.losses$cas = pmax(tips.ratite.losses$casCas,tips.ratite.losses$droNov)
tips.ratite.losses$clades = tips.ratite.losses$anoDid + tips.ratite.losses$strCam +  tips.ratite.losses$rhea + tips.ratite.losses$kiwi +  tips.ratite.losses$cas
tips.ratite.losses$total = rowSums(tips.ratite.losses[,c(1:9)])

#add ce.info
ratites.conv.losses.final.bayesian<-merge(tips.ratite.losses, ce.info[,c(2,4,5,7,8)], by.x=0, by.y=1)
write.table(ratites.conv.losses.final.bayesian, file="final_bayesian_conv_Jan2017.txt", quote=F, row.names=F, sep="\t")

#ages - note that ratite convergent seem to be older than all CNEEs

nodekey=data.frame(node=factor(c("vertebrate","node_112","node_108","node_107","node_102","node_99","node_95","node_91","node_79","node_59","node_58"), levels=c("vertebrate","node_112","node_108","node_107","node_102","node_99","node_95","node_91","node_79","node_59","node_58")), name=c("vertebrate", "tetrapods","amniotes","reptiles","turtles","Archosaurs","Birds","Neognaths","Galloanseriformes","Galliformes","Chicken"))


ratites.conv.losses.final.bayesian=merge(nodekey,ratites.conv.losses.final.bayesian, by.y="ancestralnode", by.x="node")
ce.info<-merge(nodekey, ce.info, by.y="ancestralnode", by.x="node")

ratite.origin.enrich<-c(table(ratites.conv.losses.final.bayesian$node[ratites.conv.losses.final.bayesian$clades>1])/length(ratites.conv.losses.final.bayesian$node[ratites.conv.losses.final.bayesian$clades>1]))[1:7]/c(table(ce.info$node[ce.info$node %in% c("node_102", "node_107", "node_108", "node_112", "node_95", 'node_99', "vertebrate")])/length(ce.info$node[ce.info$node %in% c("node_102", "node_107", "node_108", "node_112", "node_95", 'node_99', "vertebrate")]))[1:7]
plot(log2(ratite.origin.enrich), type="b", col="blue", lwd=2, pch=16, ylim=c(-3,3), xaxt="n", bty="n", xlab="", ylab="Relative RAR origination rate", las=2)
axis(side=1, at=c(1,2,3,4,5,6,7), labels=nodekey$name[1:7], las=2)
