library(ape)

setwd("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/blastp")
blastp<-read.table("all_besthit.out", header=F)
names(blastp)=c("file", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart","send", "evalue", "bitscore")

#add species column
blastp$species = sub("_galGal_blastp.out", "", blastp$file, fixed=T)
blastp$spclass = "Neognath"
blastp$spclass[blastp$species %in% c("eudEle", "notPer", "tinGut", "cryCin")] = "Tinamou"
blastp$spclass[blastp$species %in% c("aptRow", "aptOwe", "aptHaa", "droNov", "casCas", "rhePen", "rheAme", "strCam")] = "Ratite"
blastp$spclass[blastp$species %in% c("allMis", "anoCar", "chrPic")] = "Reptile"

#protein species
protsp<-unique(blastp$species)

#compute presence/absence matrix
prot_len = blastp[blastp$species=="galGal", c("qseqid", "qend")]
names(prot_len)[2]="qlen"
blastp<-merge(blastp, prot_len, by="qseqid")

blastp$pres_strict = ifelse(blastp$pident >= 50 & blastp$length/blastp$qlen > 0.70, 1, 0)

prot_pres_strict = data.frame(prot=blastp$qseqid[blastp$species=="galGal" & blastp$pres_strict==1])

for (sp in protsp) {
  prot_pres_strict[,sp] = as.numeric(prot_pres_strict$prot %in% blastp[blastp$pres_strict == 1 & blastp$species == sp, c("qseqid")])
  }
#get tree
sptree<-read.tree("../final_neut_tree.nwk")

#remove tips not in protsp
tips_to_drop<-sptree$tip.label[!(sptree$tip.label %in% protsp)]
prot_tree<-drop.tip(sptree, tips_to_drop)
plot(prot_tree)

#clean up to require presence at base of birds, and to remove outgroups
prot_pres_good = apply(prot_pres_strict[,c("eudEle", "notPer", "tinGut", "cryCin","aptRow", "aptOwe", "aptHaa", "droNov", "casCas", "rhePen", "rheAme", "strCam")], 1, max)
prot_pres_use = subset(prot_pres_strict, prot_pres_good == 1, select=c(-anoCar, -allMis, -chrPic))

compute_losses <- function(pres_line, orig_tree) {
  if (sum(pres_line)==length(pres_line)) { return(0) }
  newtree<-drop.tip(orig_tree, names(pres_line[pres_line==0]), subtree=T)
  neognath<-"taeGut-galGal"
  palaeo<-"aptHaa-strCam"
  tinamou<-"cryCin-eudEle"
  
  pal_loss<-tryCatch(sum(grepl("_tip", extract.clade(newtree, palaeo)$tip.label)),error = function(c) return(1))
  tin_loss<-tryCatch(sum(grepl("_tip", extract.clade(newtree, tinamou)$tip.label)),error = function(c) return(1))
  rat_loss = pal_loss - tin_loss
  return(rat_loss)
}

bootstrap_loss <- function(pres_line, orig_tree, nreps=100) {
  real_ct <- compute_losses(pres_line, orig_tree)
  perm_ct<-numeric(nreps)
  for (rep in 1:nreps) {
    shuffled_line<-pres_line
    names(shuffled_line)=sample(names(pres_line))
    perm_ct[rep]<-compute_losses(shuffled_line,orig_tree)
  }
  list(realct = real_ct, perm_ct = perm_ct, lowerp=sum(perm_ct >= real_ct)/nreps, upperp=sum(perm_ct <= real_ct)/nreps)
}

##EDIT BELOW HERE##

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
