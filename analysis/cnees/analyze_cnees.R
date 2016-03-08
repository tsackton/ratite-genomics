#prep for permutations and get 'real' results

#read in previous outputs
cnee.long<-read.table("cnee_long_info.tsv", sep="\t", header=T, stringsAsFactors=F)
ce.annot <- read.table("./final_beds/ce_annotation.tsv", header=T, stringsAsFactors = F)

#make total accelerated narrow and total accelerated broad cts
cnee.long$total.accel.br =apply(cnee.long[,c("aptHaa", "aptOwe", "aptRow", "casCas", "droNov", "strCam", "rheAme", "rhePen", "Rhea", "Casuar", "Kiwi", "ratite")], 1, function(x) sum(x <= 0.05))
cnee.long$total.accel.na = apply(cnee.long[,c("aptHaa", "aptOwe", "aptRow", "casCas", "droNov", "strCam", "rheAme", "rhePen", "Rhea", "Casuar", "Kiwi", "ratite")], 1, function(x) sum(x <= 0.01))

#make column for accelerated in at least one ratite
cnee.long$accel.broad = 0
cnee.long$accel.broad[cnee.long$total.accel.br > 0 & cnee.long$tinPres >= 2 & cnee.long$tinamou > 0.1 & cnee.long$AvesPres >= 23] = 1
cnee.long$accel.narrow = 0
cnee.long$accel.narrow[cnee.long$total.accel.na > 0 & cnee.long$tinamou > 0.25 & cnee.long$AvesPres == 35] = 1

#make column for each clade
cnee.long$casuar.br = 0
cnee.long$casuar.br[cnee.long$Casuar <= 0.05 & cnee.long$tinPres >= 2 & cnee.long$tinamou > 0.1 & cnee.long$AvesPres >= 23] = 1

cnee.long$strcam.br = 0
cnee.long$strcam.br[cnee.long$strCam <= 0.05 & cnee.long$tinPres >= 2 & cnee.long$tinamou > 0.1 & cnee.long$AvesPres >= 23] = 1

cnee.long$kiwi.br = 0
cnee.long$kiwi.br[cnee.long$Kiwi <= 0.05 & cnee.long$tinPres >= 2 & cnee.long$tinamou > 0.1 & cnee.long$AvesPres >= 23] = 1

cnee.long$rhea.br = 0
cnee.long$rhea.br[cnee.long$Rhea <= 0.05 & cnee.long$tinPres >= 2 & cnee.long$tinamou > 0.1 & cnee.long$AvesPres >= 23] = 1

#also do strict column for each clade
cnee.long$casuar.na = 0
cnee.long$casuar.na[cnee.long$Casuar <= 0.01  & cnee.long$tinamou > 0.25 & cnee.long$AvesPres == 35] = 1

cnee.long$strcam.na = 0
cnee.long$strcam.na[cnee.long$strCam <= 0.01  & cnee.long$tinamou > 0.25 & cnee.long$AvesPres == 35] = 1

cnee.long$kiwi.na = 0
cnee.long$kiwi.na[cnee.long$Kiwi <= 0.01  & cnee.long$tinamou > 0.25 & cnee.long$AvesPres == 35] = 1

cnee.long$rhea.na = 0
cnee.long$rhea.na[cnee.long$Rhea <= 0.01  & cnee.long$tinamou > 0.25 & cnee.long$AvesPres == 35] = 1

cnee.long$clade.ct.br = apply(cnee.long[,c("rhea.br", "kiwi.br", "strcam.br", "casuar.br")], 1, sum)
cnee.long$clade.ct.na = apply(cnee.long[,c("rhea.na", "kiwi.na", "strcam.na", "casuar.na")], 1, sum)


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


