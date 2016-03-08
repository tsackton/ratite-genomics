library(plyr)

#prep for permutations and get 'real' results

#read in ce annotation
ce.annot <- read.table("../phast/cnee_long_info.tsv", header=T, stringsAsFactors = F)
#cnee.data<-

#now want to test robustness to tree model and analysis type, so need four subsets -- original, tree 1, tree 2, tree 3
#data comes from add_acceleration_moa.R
#need the moa.wide data frame but everything else is created

moa.for.robustness<-moa.wide[,c("name", "model", "anoDid", "strCam", "Casuar", "Kiwi", "Rhea", "tinamou", "aptHaa", "aptOwe", "aptRow", "casCas", "droNov", "rheAme", "rhePen", "allRatite", "basalPaleo")]
moa.orig<-ce.annot[,c("id", "ratite", "strCam", "Casuar", "Kiwi", "Rhea", "tinamou", "aptHaa", "aptOwe", "aptRow", "casCas", "droNov", "rheAme", "rhePen")]
names(moa.orig)[1]="name"
names(moa.orig)[2]="allRatite"
moa.orig$anoDid=NA
moa.orig$basalPaleo=NA
moa.orig$model="original"
moa.robustness<-rbind(moa.for.robustness, moa.orig)
moa.robustness<-merge(moa.robustness, cnee.data[,c("id", "best_ens", "best_ncbi", "tinPres", "AvesPres", "Freq")], by.y="id", by.x="name") 

#broad column on moa robustness data
moa.robustness$casuar.br = 0
moa.robustness$casuar.br[moa.robustness$Casuar <= 0.05  & moa.robustness$tinamou > 0.1 & moa.robustness$tinPres >= 2  & moa.robustness$AvesPres >= 23 & (moa.robustness$basalPaleo > 0.1 | is.na(moa.robustness$basalPaleo)) & moa.robustness$Freq < 5] = 1

moa.robustness$strcam.br = 0
moa.robustness$strcam.br[moa.robustness$strCam <= 0.05 & moa.robustness$tinamou > 0.1 & moa.robustness$tinPres >= 2 & moa.robustness$AvesPres >= 23 & (moa.robustness$basalPaleo > 0.1 | is.na(moa.robustness$basalPaleo)) & moa.robustness$Freq < 5] = 1

moa.robustness$kiwi.br = 0
moa.robustness$kiwi.br[moa.robustness$Kiwi <= 0.05 & moa.robustness$tinamou > 0.1 & moa.robustness$tinPres >= 2 & moa.robustness$AvesPres >= 23 & (moa.robustness$basalPaleo > 0.1 | is.na(moa.robustness$basalPaleo)) & moa.robustness$Freq < 5] = 1

moa.robustness$rhea.br = 0
moa.robustness$rhea.br[moa.robustness$Rhea <= 0.05  & moa.robustness$tinamou > 0.1 & moa.robustness$tinPres >= 2 & moa.robustness$AvesPres >= 23 & (moa.robustness$basalPaleo > 0.1 | is.na(moa.robustness$basalPaleo)) & moa.robustness$Freq < 5] = 1

moa.robustness$moa.br = 0
moa.robustness$moa.br[!is.na(moa.robustness$anoDid) & moa.robustness$anoDid <= 0.05 & moa.robustness$tinamou > 0.1 & moa.robustness$tinPres >= 2 & moa.robustness$AvesPres >= 23 & (moa.robustness$basalPaleo > 0.1 | is.na(moa.robustness$basalPaleo)) & moa.robustness$Freq < 5] = 1

moa.robustness$br.ct = moa.robustness$casuar.br + moa.robustness$strcam.br + moa.robustness$kiwi.br + moa.robustness$rhea.br  + moa.robustness$moa.br

#narrow column on moa robustness data
moa.robustness$casuar.na = 0
moa.robustness$casuar.na[moa.robustness$Casuar <= 0.01  & moa.robustness$tinamou > 0.25 & moa.robustness$AvesPres >= 35 & (moa.robustness$basalPaleo > 0.25 | is.na(moa.robustness$basalPaleo)) & moa.robustness$Freq < 2] = 1

moa.robustness$strcam.na = 0
moa.robustness$strcam.na[moa.robustness$strCam <= 0.01 & moa.robustness$tinamou > 0.25 & moa.robustness$AvesPres >= 35 & (moa.robustness$basalPaleo > 0.25 | is.na(moa.robustness$basalPaleo)) & moa.robustness$Freq < 2] = 1

moa.robustness$kiwi.na = 0
moa.robustness$kiwi.na[moa.robustness$Kiwi <= 0.01 & moa.robustness$tinamou > 0.25 & moa.robustness$AvesPres >= 35 & (moa.robustness$basalPaleo > 0.25 | is.na(moa.robustness$basalPaleo)) & moa.robustness$Freq < 2] = 1

moa.robustness$rhea.na = 0
moa.robustness$rhea.na[moa.robustness$Rhea <= 0.01  & moa.robustness$tinamou > 0.25 & moa.robustness$AvesPres >= 35 & (moa.robustness$basalPaleo > 0.25 | is.na(moa.robustness$basalPaleo)) & moa.robustness$Freq < 2] = 1

moa.robustness$moa.na = 0
moa.robustness$moa.na[!is.na(moa.robustness$anoDid) & moa.robustness$anoDid <= 0.01 & moa.robustness$tinamou > 0.25 & moa.robustness$AvesPres >= 35 & (moa.robustness$basalPaleo > 0.25 | is.na(moa.robustness$basalPaleo)) & moa.robustness$Freq < 2] = 1

moa.robustness$na.ct = moa.robustness$casuar.na + moa.robustness$strcam.na + moa.robustness$kiwi.na + moa.robustness$rhea.na  + moa.robustness$moa.br

#make total accelerated narrow and total accelerated broad cts
moa.robustness$total.accel.br =apply(moa.robustness[,c("aptHaa", "aptOwe", "aptRow", "casCas", "droNov", "strCam", "rheAme", "rhePen", "Rhea", "Casuar", "Kiwi", "allRatite", "anoDid")], 1, function(x) sum(x <= 0.05, na.rm=T))
moa.robustness$total.accel.na = apply(moa.robustness[,c("aptHaa", "aptOwe", "aptRow", "casCas", "droNov", "strCam", "rheAme", "rhePen", "Rhea", "Casuar", "Kiwi", "allRatite", "anoDid")], 1, function(x) sum(x <= 0.01, na.rm=T))

#make column for accelerated in at least one ratite
moa.robustness$accel.broad = 0
moa.robustness$accel.broad[moa.robustness$total.accel.br > 0 & moa.robustness$tinPres >= 2 & moa.robustness$tinamou > 0.1 & moa.robustness$AvesPres >= 23 & (moa.robustness$basalPaleo > 0.1 | is.na(moa.robustness$basalPaleo)) & moa.robustness$Freq < 5] = 1
moa.robustness$accel.narrow = 0
moa.robustness$accel.narrow[moa.robustness$total.accel.na > 0 & moa.robustness$tinamou > 0.25 & moa.robustness$AvesPres == 35 & (moa.robustness$basalPaleo > 0.25 | is.na(moa.robustness$basalPaleo)) & moa.robustness$Freq < 2] = 1

#convert to something useful, based on broad results, for Phil to analyze
moa.conv.results<-spread(moa.robustness[,c("name", "model", "best_ens", "best_ncbi", "br.ct")], model, br.ct)

#generate sum of moa accelerations for moa final dataset
moa.accel=ddply(moa.robustness[,c("name", "moa.br")], .(name), summarize, moa.ct=sum(moa.br)) 
moa.conv.results<-merge(moa.conv.results, moa.accel)

write.table(moa.conv.results, file="final_convergence_results.tsv", sep="\t", quote=F, col.names=T, row.names=F)
write.table(moa.robustness, file="raw_convergence_results.tsv", sep="\t", quote=F, col.names=T, row.names=F)

#make a good cnee subset for remaining analysis and for permutations
cnee.sub<-subset(moa.robustness, Freq < 5 & AvesPres >= 23)

##refactor below##

#get real counts for each gene, using NCBI_best or ens_best
ens.genes.broad = table(cnee.long$best_ens, cnee.long$accel.broad)
ens.genes.broad = as.data.frame.matrix(ens.genes.broad)
names(ens.genes.broad) = c("cons", "accel")
ens.genes.broad$total = ens.genes.broad$cons + ens.genes.broad$accel
ens.genes.broad = subset(ens.genes.broad, rownames(ens.genes.broad)!="none")

ens.genes.narrow = table(cnee.long$best_ens, cnee.long$accel.narrow)
ens.genes.narrow = as.data.frame.matrix(ens.genes.narrow)
names(ens.genes.narrow) = c("cons", "accel")
ens.genes.narrow$total = ens.genes.narrow$cons + ens.genes.narrow$accel
ens.genes.narrow = subset(ens.genes.narrow, rownames(ens.genes.narrow)!="none")

ncbi.genes.broad = table(cnee.long$best_ncbi, cnee.long$accel.broad)
ncbi.genes.broad = as.data.frame.matrix(ncbi.genes.broad)
names(ncbi.genes.broad) = c("cons", "accel")
ncbi.genes.broad$total = ncbi.genes.broad$cons + ncbi.genes.broad$accel
ncbi.genes.broad = subset(ncbi.genes.broad, rownames(ncbi.genes.broad)!="none")

ncbi.genes.narrow = table(cnee.long$best_ncbi, cnee.long$accel.narrow)
ncbi.genes.narrow = as.data.frame.matrix(ncbi.genes.narrow)
names(ncbi.genes.narrow) = c("cons", "accel")
ncbi.genes.narrow$total = ncbi.genes.narrow$cons + ncbi.genes.narrow$accel
ncbi.genes.narrow = subset(ncbi.genes.narrow, rownames(ncbi.genes.narrow)!="none")

#now write out output for permutations
cnee.for.perm<-cnee.long[,c(13,14,77:86)]

write.table(cnee.for.perm, file="cnee.for.perm", sep="\t", quote=F, row.names=F)


