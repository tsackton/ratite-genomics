#analyze cons_scale results

setwd("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/analysis")
csrows<-scan("cons_scale_ver2_rows", what="character")
cscols<-scan("cons_scale_ver2_header", what="character")
cons<-read.table("cons_scale_ver2_elem_Z.txt", stringsAsFactors=F, header=F, row.names=csrows, col.names=cscols)
head(cons)
ratite<-c("aptRow", "aptOwe", "aptHaa", "casCas", "droNov", "rheAme", "rhePen", "anoDid", "strCam", "aptHaa.aptOwe", "aptHaa.aptRow", "casCas.droNov", "rheAme.rhePen")
nonbird<-c("allMis", "allSin", "croPor", "gavGan", "chrPic", "cheMyd", "anoCar", "taeGut.aptHaa", "allMis.allSin", "croPor.gavGan", "allMis.croPor", "taeGut.allMis", "chrPic.cheMyd", "taeGut.chrPic", "taeGut.anoCar")
paleo<-c("cryCin", "tinGut", "eudEle", "notPer", "aptHaa.casCas", "aptHaa.rheAme", "cryCin.tinGut", "eudEle.notPer", "cryCin.eudEle", "cryCin.anoDid", "aptHaa.cryCin", "aptHaa.strCam")
penguin<-c("pygAde", "aptFor", "aptFor.pygAde")
neognath<-colnames(cons)[!(colnames(cons) %in% nonbird | colnames(cons) %in% ratite | colnames(cons) %in% paleo | colnames(cons) %in% penguin)]
cons$ratite<-apply(cons[,ratite], 1, function(x) sum(x != 1))
cons$nonbird<-apply(cons[,nonbird], 1, function(x) sum(x != 1))
cons$penguin<-apply(cons[,penguin], 1, function(x) sum(x != 1))
cons$paleo<-apply(cons[,paleo], 1, function(x) sum(x != 1))
cons$neognath<-apply(cons[,neognath], 1, function(x) sum(x != 1))
birdonly<-subset(cons, cons$taeGut.aptHaa == 1)

random<-c("taeGut","nipNip", "falPer", "colLiv")
random.set<-cons[,random]
random.set$ct = apply(random.set, 1, function(x) sum(x != 1))
ratite.sub<-c("aptRow", "casCas", "rheAme", "strCam")
ratite.set<-cons[,ratite.sub]
ratite.set$ct = apply(ratite.set, 1, function(x) sum(x != 1))

highcons<-subset(birdonly, paleo==0 & neognath < 2)
highcons2<-subset(birdonly, paleo+ratite < 2)

#temp
cands<-read.table("cnee_permutations/ratite_candidates_final.tsv", header=T, sep="\t")
cands<-merge(cands, cons, all.x=T, all.y=F, by.x="name", by.y="row.names")
cands2<-subset(cands, neognath < 2 & paleo==0)

barplot(table(highcons2$neognath)/length(highcons$neognath), names=c(round(seq(0,1,by=1/42),2)), las=2)

barplot(table(highcons$ratite)/length(highcons$ratite), las=1, xlab="# of ratite lineages where Z != conserved")
barplot(table(cands$ratite)/length(cands$ratite), las=1, xlab="# of ratite lineages where Z != conserved")

