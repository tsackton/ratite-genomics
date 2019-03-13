# read in annotation
anno <- read.table("cnees.galgal4.annotation", stringsAsFactors = F)
anno = anno[grepl("(DACH1|NF1B|NPAS3)",anno[,2]),]


# read in GC
# GC<-fread("newData/allspecies_cnee_concat_GC.txt", stringsAsFactors = F)
# write.table(t(GC), "newData/allspecies_cnee_concat_GC2.txt",sep="\t", quote=F, col.names = F)

Tot<-fread("GC/allspecies_cnee_concat_Tot2.txt")
GC<-fread("GC/allspecies_cnee_concat_GC2.txt")


# read in species name
species_name <- read.table("species_names.txt", stringsAsFactors = F)


# read in accelerated elements
score <- read.table("Gain/Combined_elem_lik.txt", sep="\t", header=T,stringsAsFactors = F)
postZ <- fread("Gain/Combined_post_Z_M2.txt", header=F)

postZ_acc <- postZ[,.SD, .SDcols=c(1,seq(10, 178, by = 4))] # get ncol(postZ), Prob(Z==2)
setkey(postZ_acc,V1)
colnames(postZ_acc)[2:44] <- species_name[1:43,1]

postZ_cons <- postZ[,.SD, .SDcols=c(1,seq(9, 178, by = 4))] # get ncol(postZ), Prob(Z==2)
setkey(postZ_cons,V1)
colnames(postZ_cons)[2:44] <- species_name[1:43,1]


#### Compare GC content for accelerated vs. conserved branches ####

# get accelerated elements
sel <- score[which(score$loglik_Full - score$loglik_Null > 10),"ID"]
# get accelerated elements within HAR
sel1 <- intersect(anno[,1], sel) #100
# get RAR
sel2 <- score[score$logBF1 > 10 & score$logBF2 > 1, "ID"]

# reorder columns of GC to be the same as postZ (only birds)
GC1 = GC[, c("Elements", species_name[1:36,1]),with=F] 
Tot1 = Tot[, c("Elements",species_name[1:36,1]),with=F]
setkey(GC1,Elements)
setkey(Tot1,Elements)

GC_sel = rowSums(GC1[sel,-1])/rowSums(Tot1[sel,-1])
GC_sel1 = rowSums(GC1[sel1,-1])/rowSums(Tot1[sel1,-1])
GC_sel2 = rowSums(GC1[sel2,-1])/rowSums(Tot1[sel2,-1])

#get accelerated species
selPostZ = (postZ_acc[sel,2:37, with=F] > 0.7) #44
GC_sel_acc = rowSums(GC1[sel,-1] * selPostZ)/rowSums(Tot1[sel,-1] * selPostZ)

# get conserved species
selPostZ = (postZ_cons[sel,2:37, with=F] > 0.7)
GC_sel_cons = rowSums(GC1[sel,-1] * selPostZ)/rowSums(Tot1[sel,-1] * selPostZ)

#get accelerated species for HAR
selPostZ = (postZ_acc[sel1,2:37, with=F] > 0.7)
GC_sel1_acc = rowSums(GC1[sel1,-1] * selPostZ)/rowSums(Tot1[sel1,-1] * selPostZ)

# get conserved species for HAR
selPostZ = (postZ_cons[sel1,2:37, with=F] > 0.7)
GC_sel1_cons = rowSums(GC1[sel1,-1] * selPostZ)/rowSums(Tot1[sel1,-1] * selPostZ)


#get accelerated species for RAR
selPostZ = (postZ_acc[sel2,2:37, with=F] > 0.7)
GC_sel2_acc = rowSums(GC1[sel2,-1] * selPostZ)/rowSums(Tot1[sel2,-1] * selPostZ)

# get conserved species for RAR
selPostZ = (postZ_cons[sel2,2:37, with=F] > 0.7)
GC_sel2_cons = rowSums(GC1[sel2,-1] * selPostZ)/rowSums(Tot1[sel2,-1] * selPostZ)



#x1 = c("32918 CNEEs", "111 HAR CNEEs", "2759 RAR CNEEs")
x1 = c("accelerated", "HAR", "RAR")
x2 = c("All Birds",  "Conserved bird lineages","Accelerated bird lineages")
label = paste(rep(x2,3),rep(x1, each=3), sep="\n")

pdf("Figures/GC/GC_HAR_RAR_boxplot_acceleration_onlybirds2.pdf", width=18)
#par(cex=1.2)
boxplot(list(GC_sel, GC_sel_cons, GC_sel_acc,GC_sel1, GC_sel1_cons, GC_sel1_acc,GC_sel2, GC_sel2_cons, GC_sel2_acc),
        names= label , outline = F, col = rep(c(3,4,7), each=3),notch=T, ylab = "GC content (%GC)", main="accelerated CNEEs")
abline(h=sum(GC1[,-1])/sum(Tot1[,-1]) , col=2, lty=2 ) #0.3776348
abline(h=sum(GC1[anno[,1],-1])/sum(Tot1[anno[,1],-1]), col=5, lty=2 ) #  0.3529689
legend("topleft", c("all CNEEs", "HAR CNEEs"),title="mean GC", lty=2, col=c(2,5))
dev.off()



