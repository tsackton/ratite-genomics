setwd("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/blastp")
blastp<-read.table("all_besthit.out", header=F)
names(blastp)=c("file", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart","send", "evalue", "bitscore")

#add species column
blastp$species = sub("_galGal_blastp.out", "", blastp$file, fixed=T)
blastp$spclass = "Neognath"
blastp$spclass[blastp$species %in% c("eudEle", "notPer", "tinGut", "cryCin")] = "Tinamou"
blastp$spclass[blastp$species %in% c("aptRow", "aptOwe", "aptHaa", "droNov", "casCas", "rhePen", "rheAme", "strCam")] = "Ratite"
blastp$spclass[blastp$species %in% c("allMis", "anoCar", "chrPic")] = "Reptile"

#no filtering, count of species present
spcount<-as.data.frame(table(blastp$qseqid))
table(spcount$Freq)

#some percent ID > 70%
spcount.70<-as.data.frame(table(blastp$qseqid[blastp$pident >= 70]))
table(spcount.70$Freq)

#species table
spmat<-as.data.frame(table(blastp$qseqid, blastp$species))
spclassmat<-as.data.frame(table(blastp$qseqid, blastp$spclass))

#unstack
spclass.wide<-reshape(spclassmat, timevar="Var2", idvar="Var1", direction="wide")
table(spclass.wide$Freq.Ratite, spclass.wide$Freq.Tinamou)

spclass.wide$ratite = spclass.wide$Freq.Ratite / 8
spclass.wide$tinamou = spclass.wide$Freq.Tinamou / 4
spclass.wide$neognath = spclass.wide$Freq.Neognath / 27

plot(spclass.wide$ratite ~ spclass.wide$neognath)

#tinamou present
spclass.wide.tp<-subset(spclass.wide, tinamou==1)
table(spclass.wide.tp$Freq.Ratite, spclass.wide.tp$Freq.Neognath)

spclass.wide.tp$ratite.p = NA
spclass.wide.tp$ratite.or = NA

for (i in 1:length(spclass.wide.tp$Freq.Neognath)) {
  spclass.wide.tp$ratite.p[i] = fisher.test(matrix(c(spclass.wide.tp$Freq.Ratite[i], 8-spclass.wide.tp$Freq.Ratite[i], spclass.wide.tp$Freq.Neognath[i], 27-spclass.wide.tp$Freq.Neognath[i]),nrow=2))$p.val
  spclass.wide.tp$ratite.or[i] = fisher.test(matrix(c(spclass.wide.tp$Freq.Ratite[i], 8-spclass.wide.tp$Freq.Ratite[i], spclass.wide.tp$Freq.Neognath[i], 27-spclass.wide.tp$Freq.Neognath[i]),nrow=2))$est
}

spclass.wide.tp$qval = p.adjust(spclass.wide.tp$ratite.p, method="fdr")
spclass.wide.tp$bon = p.adjust(spclass.wide.tp$ratite.p, method="holm")

spclass.wide.tp[spclass.wide.tp$ratite.p < 0.01 & spclass.wide.tp$ratite.or < 1,]
ratite.loss<-spclass.wide.tp$Var1[spclass.wide.tp$ratite.p < 0.01 & spclass.wide.tp$ratite.or < 1]
ratite.loss

blastp[blastp$qseqid=="ENSGALP00000017755.1",]
blastp[blastp$qseqid=="ENSGALP00000011548.3",]
blastp[blastp$qseqid=="ENSGALP00000013094.2",]
