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
neognath<-colnames(cons)[!(colnames(cons) %in% nonbird | colnames(cons) %in% ratite | colnames(cons) %in% paleo | colnames(cons) %in% penguin)]
allbirds<-c(neognath,penguin,ratite,paleo)
cons$ratite<-apply(cons[,ratite], 1, function(x) sum(x != 1))
cons$nonbird<-apply(cons[,nonbird], 1, function(x) sum(x != 1))
cons$penguin<-apply(cons[,penguin], 1, function(x) sum(x != 1))
cons$paleo<-apply(cons[,paleo], 1, function(x) sum(x != 1))
cons$neognath<-apply(cons[,neognath], 1, function(x) sum(x != 1))
#birdonly<-subset(cons, cons$taeGut.aptHaa == 1)

#random<-c("taeGut","nipNip", "falPer", "colLiv")
#random.set<-cons[,random]
#random.set$ct = apply(random.set, 1, function(x) sum(x != 1))
#ratite.sub<-c("aptRow", "casCas", "rheAme", "strCam")
#ratite.set<-cons[,ratite.sub]
#ratite.set$ct = apply(ratite.set, 1, function(x) sum(x != 1))

#highcons<-subset(birdonly, paleo==0 & neognath < 2)
#highcons2<-subset(birdonly, paleo+ratite < 2)

#temp
#cands<-read.table("cnee_permutations/ratite_candidates_final.tsv", header=T, sep="\t")
#cands<-merge(cands, cons, all.x=T, all.y=F, by.x="name", by.y="row.names")
#cands2<-subset(cands, neognath < 2 & paleo==0)

#barplot(table(highcons2$neognath)/length(highcons$neognath), names=c(round(seq(0,1,by=1/42),2)), las=2)

#barplot(table(highcons$ratite)/length(highcons$ratite), las=1, xlab="# of ratite lineages where Z != conserved")
#barplot(table(cands$ratite)/length(cands$ratite), las=1, xlab="# of ratite lineages where Z != conserved")


cons$ratite.test<-apply(cons[,c("aptRow", "casCas", "rheAme", "strCam","anoDid")], 1, function(x) sum(x != 1))
cons$ratite.3.test<-apply(cons[,c("aptRow","strCam","anoDid")], 1, function(x) sum(x != 1))
cons$ratite.control<-apply(cons[,allbirds[!allbirds %in% ratite]], 1, function(x) sum(x != 1))
cons$random1.test<-apply(cons[,c("colLiv", "chaVoc", "halLeu", "taeGut", "nipNip")], 1, function(x) sum(x != 1))
cons$random1.control<-apply(cons[,allbirds[!allbirds %in% c("colLiv", "chaVoc", "halLeu", "taeGut", "nipNip")]], 1, function(x) sum(x != 1))
cons$vl.test<-apply(cons[,c("calAnn", "melUnd", "ficAlb")], 1, function(x) sum(x != 1))
cons$vl.control<-apply(cons[,allbirds[!allbirds %in% c('calAnn', 'corBra', 'serCan', 'geoFor', 'melUnd', 'pseHum', 'taeGut', 'ficAlb')]], 1, function(x) sum(x != 1))
cons$random2.test<-apply(cons[,c("falPer", "picPub", "lepDis")], 1, function(x) sum(x != 1))
cons$random2.control<-apply(cons[,allbirds[!allbirds %in% c("falPer", "picPub", "lepDis")]], 1, function(x) sum(x != 1))

barplot(t(matrix(c(table(cons$ratite.test[cons$ratite.control==0 & cons$ratite.test>1]), c(table(cons$random1.test[cons$random1.control==0 & cons$random1.test>1]),0,0)),ncol=2)),beside=T,col=c("red", "blue"),names=c("2", "3", "4", "5"), ylab="# of Elements", xlab="# of Lineages", main="Elements Uniquely Accelerated in Target Lineages")
legend("topright", legend=c("Ratite", "Random Control"), col=c("red", "blue"), pch=15, bty="n")


barplot(t(matrix(c(table(cons$vl.test[cons$vl.control==0 & cons$vl.test>1]), c(table(cons$random2.test[cons$random2.control==0 & cons$random2.test>1]),0)),ncol=2)),beside=T,col=c("darkgreen", "blue"),names=c("2", "3"), ylab="# of Elements", xlab="# of Lineages", main="Elements Uniquely Accelerated in Target Lineages")
legend("topright", legend=c("Vocal Learners", "Random Control"), col=c("darkgreen", "blue"), pch=15, bty="n")

barplot(t(matrix(c(table(cons$ratite.3.test[cons$ratite.control==0 & cons$ratite.3.test>1]), c(table(cons$random2.test[cons$random2.control==0 & cons$random2.test>1]),0)),ncol=2)),beside=T,col=c("red", "blue"),names=c("2", "3"), ylab="# of Elements", xlab="# of Lineages", main="Elements Uniquely Accelerated in Target Lineages")
legend("topright", legend=c("Ratite", "Random Control"), col=c("red", "blue"), pch=15, bty="n")

barplot(t(matrix(c(table(cons$ratite.3.test[cons$ratite.control==0 & cons$ratite.3.test>1]), table(cons$vl.test[cons$vl.control==0 & cons$vl.test>1]), c(table(cons$random2.test[cons$random2.control==0 & cons$random2.test>1]),0)),ncol=3)),beside=T,col=c("red", "darkgreen", "blue"),names=c("2", "3"), ylab="# of Elements", xlab="# of Lineages", main="Elements Uniquely Accelerated in Target Lineages")
legend("topright", legend=c("Ratite", "Vocal Learners", "Random Control"), col=c("red", "darkgreen", "blue"), pch=15, bty="n")


barplot(t(matrix(c(table(cons$ratite.3.test[cons$ratite.control<3 & cons$ratite.3.test>1]), table(cons$vl.test[cons$vl.control<3 & cons$vl.test>1]), table(cons$random2.test[cons$random2.control<3 & cons$random2.test>1])),ncol=3)),beside=T,col=c("red", "darkgreen", "blue"),names=c("2", "3"), ylab="# of Elements", xlab="# of Lineages", main="Elements Uniquely Accelerated in Target Lineages")
legend("topright", legend=c("Ratite", "Vocal Learners", "Random Control"), col=c("red", "darkgreen", "blue"), pch=15, bty="n")

barplot(t(matrix(c(table(cons$ratite.test[cons$ratite.control<3 & cons$ratite.test>1]), c(table(cons$random1.test[cons$random1.control<3 & cons$random1.test>1]),0)),ncol=2)),beside=T,col=c("red", "blue"),names=c("2", "3", "4", "5"), ylab="# of Elements", xlab="# of Lineages", main="Elements Uniquely Accelerated in Target Lineages")
legend("topright", legend=c("Ratite", "Random Control"), col=c("red", "blue"), pch=15, bty="n")

cons$ratite.conv<-0
cons$ratite.conv[cons$ratite>2 & cons$ratite.control < 2]<-1

cons$vl.conv<-0
cons$vl.conv[cons$vl.control < 2 & cons$vl.test>1]<-1

cons$origin.crude=NA
cons$origin.crude[cons$taeGut.aptHaa==0]="in.birds"
cons$origin.crude[cons$taeGut.aptHaa==1]="base.birds"
cons$origin.crude[cons$taeGut.allMis==1]="base.arch"
cons$origin.crude[cons$taeGut.chrPic==1]="rept"

prop.table(table(cons$ratite.conv[cons$origin.crude != "in.birds"], cons$origin.crude[cons$origin.crude != "in.birds"]),1)

acc.elements<-rownames(cons[cons$ratite.conv==1,])
acc.elem.vl<-rownames(cons[cons$vl.conv==1,])

#compare to previous
orig.cand<-read.table("../cnee_permutations/ratite_candidates_final.tsv", header=T, sep="\t")
orig.cand$cons.scale = orig.cand$name %in% acc.elements
#this doesn't work because I need to read NCBI key in
#ce.info<-read.csv("~/Dropbox/Work/EggCNEEs/October 2016/gg.final.annotated.csv")
#ce.info$ratite.conv = ce.info$gg.cnee.id %in% acc.elements
#ce.info<-merge(ce.info, ncbikey, by.x="best_ens", by.y="ensembl", all.x=T, all.y=F)
#ce.info$vl.conv = ce.info$gg.cnee.id %in% acc.elem.vl
#genes<-unique(ce.info$best_ens[ce.info$ratite.conv])
#background<-unique(ce.info$best_ens)
#write.table(genes, file="ratite.conv.cnee.genes", quote=F, row.names=F, col.names=F)
#write.table(ce.info$best_ens[ce.info$vl.conv], file="vl.conv.cnee.genes", quote=F, row.names=F, col.names=F)

#write.table(background, file="background.cnee.genes", quote=F, row.names=F, col.names=F)

####making final Bayesian file for Phil####
#add family level ratite data
rhea<-c("rheAme", "rhePen", "rheAme.rhePen")
kiwi<-c("aptHaa", "aptRow", "aptOwe", "aptHaa.aptOwe", "aptHaa.aptRow")
cas<-c("casCas", "droNov", "casCas.droNov")
moa<-c("anoDid")
ost<-c("strCam")

cons$rhea<-apply(cons[,rhea], 1, function(x) sum(x != 1)/3)
cons$kiwi<-apply(cons[,kiwi], 1, function(x) sum(x != 1)/5)
cons$cas<-apply(cons[,cas], 1, function(x) sum(x != 1)/3)
cons$ostrich<-apply(cons[,ost, drop=F], 1, function(x) sum(x != 1)/1)
cons$moa<-apply(cons[,moa, drop=F], 1, function(x) sum(x != 1)/1)

cons$ratite.conv<-apply(cons[,c("rhea", "kiwi", "cas", "ostrich","moa")], 1, function(x) sum(x == 1))
cons$ratite.conv[cons$paleo > 0] = 0

#subset for Phil
cons.cand<-subset(cons, select=c("ratite.conv", "kiwi", "cas", "rhea", "ostrich", "moa", "rheAme", "nonbird", "ratite", "penguin", "neognath", "paleo"))
cons.cand$rheAme[cons.cand$rheAme != "1"] = "neut/missing"
cons.cand$rheAme[cons.cand$rheAme == "1"] = "conserved"

write.table(cons.cand, file="cons_cand_2016Dec06.txt", row.names=T, sep="\t", quote=F, col.names=T)
