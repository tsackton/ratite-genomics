library(ape)
library(dplyr)
library(tidyr)
library(HyPhy) #used for gene tree species tree reconcilation to get loss numbers

setwd("/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/blastp")
blastp<-read.table("all_besthit.out", header=F, stringsAsFactors = F)
blastp <- tbl_df(blastp)
blastp <- blastp %>% 
  rename(file = V1, query_id = V2, subj_id = V3, perc_ident = V4, length = V5, mismatch = V6, gapopen = V7, query_start = V8, query_end = V9, subj_start = V10, subj_end = V11, eval = V12, bitscore = V13) %>%
  mutate(species = sub("_galGal_blastp.out", "", file, fixed=T))

#make species <-> clade key
tinsp <- c("eudEle", "notPer", "tinGut", "cryCin")
ratitesp <- c("aptRow", "aptOwe", "aptHaa", "droNov", "casCas", "rhePen", "rheAme", "strCam")
reptilesp <- c("allMis", "anoCar", "chrPic")
protsp <- blastp %>% distinct(species) %>% group_by(species) %>%
  mutate(spclass = { if(species %in% tinsp) {return("tinamou")} else if (species %in% ratitesp) {return("ratite")} else if (species %in% reptilesp) { return("reptile")} else { return ("neognath")}  })

#merge

blastp <- blastp %>% inner_join(., protsp)

#add length, also has the effect of removing a few sequences that don't blast back to chicken for some reason

blastp <- blastp %>% filter(species == "galGal") %>% 
  select(query_id, galgal_len = query_end) %>% 
  inner_join(blastp, ., by=c("query_id" = "query_id"))

#compute presence/absence matrix

blastp <- blastp %>% mutate(pres_strict = as.numeric(perc_ident >= 50 & length/galgal_len > 0.70), pres_loose = 1)

#convert to wide format

blastp_wide <- blastp %>% select(query_id, species, pres_strict) %>% spread(species, pres_strict, drop=FALSE, fill=0)

#get tree
full_species_tree<-read.tree("../final_neut_tree.nwk")

#remove tips not in protsp
tips_to_drop<-full_species_tree$tip.label[!(full_species_tree$tip.label %in% protsp$species)]
prot_tree<-drop.tip(full_species_tree, tips_to_drop)
plot(prot_tree)

#working with blastp_wide, do some cleanup

blastp_use <- blastp_wide %>% mutate(use_prot = pmax(eudEle, notPer, tinGut, cryCin, aptRow, aptOwe, aptHaa, droNov, casCas, rhePen, rheAme, strCam)) %>% filter(use_prot == 1) %>% select(-anoCar, -allMis, -chrPic) %>% gather(species, present, anaPla:tinGut)

#15083 genes inferred to be present at base of birds (i.e., present in galGal and at least 1 palaeognath)

compute_losses <- function(DF, orig_tree) {
  pres_line = unlist(DF$present)
  spec_line = unlist(DF$species[DF$present == 0])
  newtree<-drop.tip(orig_tree, spec_line)
  node_names = list("neognath" = "taeGut-galGal", "palaeo" = "aptHaa-strCam", "tinamou" = "cryCin-eudEle")
  
  #output
  output<-data.frame(clade=c("neognath", "palaeo", "tinamou"), loss=NA_integer_, br=NA_integer_)
  
  #loop over branches 
  for (i in 1:3) {
    sptree<-extract.clade(orig_tree, node_names[[i]])
    gtree<-tryCatch(extract.clade(newtree, node_names[[i]]), error=function(c) {return(0)})
    if (!is.numeric(gtree)) {
      output$loss[output$clade == names(node_names)[i]]<-recon.score(sptree, gtree)[2]
      output$br[output$clade == names(node_names)[i]]<-(length(gtree$edge.length)
    }
    else {
      output$loss[output$clade == names(node_names)[i]]<-1
      output$br[output$clade == names(node_names)[i]]<-1
    }
  }
  return(output)
}

#look at number of ratite losses

protein_loss <- blastp_use %>% group_by(query_id) %>% do(compute_losses(DF=., orig_tree=prot_tree))

protein_loss %>% select(query_id, clade, loss) %>% spread(clade, loss) %>% filter(tinamou == 0) %>% with(., table(neognath, palaeo))

ratite_treelen <- sum(extract.clade(full_species_tree, "aptHaa-strCam")$edge.length) - sum(extract.clade(full_species_tree, "cryCin-eudEle")$edge.length) 
neo_treelen <- sum(extract.clade(full_species_tree, "taeGut-galGal")$edge.length)

protein_loss %>% select(query_id, clade, loss) %>% spread(clade, loss) %>% filter(tinamou == 0, palaeo>3) %>% left_join(., blastp_wide) %>% mutate(ratite_ct = aptHaa + aptOwe + aptRow + casCas + droNov + rheAme + rhePen + strCam) %>% print.data.frame

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
