# read in annotation
anno <- read.table("newData/cnees.galgal4.annotation")
anno = anno[grepl("(DACH1|NF1B|NPAS3)",anno[,2]),]

# read in GC
GC<-fread("newData/allspecies_cnee_concat_GC2.txt")
Tot<-fread("newData/allspecies_cnee_concat_Tot2.txt")
setkey(GC,Elements)
setkey(Tot,Elements)

# read in species name
species_name <- read.table("Data/species_names3.txt", stringsAsFactors = F)
ratites <- c("strCam","rhePen","rheAme","casCas","droNov","aptRow","aptHaa","aptOwe","anoDid") 
neognathae <- species_name[1:23,1]

# read in accelerated elements
score <- read.table("Gain_1012/Combined_elem_lik.txt", sep="\t", header=T,stringsAsFactors = F)
postZ <- fread("Gain_1012/Combined_post_Z_M2.txt", header=F)

postZ_acc <- postZ[,.SD, .SDcols=c(1,seq(10, 178, by = 4))] # get ncol(postZ), Prob(Z==2)
setkey(postZ_acc,V1)
colnames(postZ_acc)[2:44] <- species_name[1:43,1]

postZ_cons <- postZ[,.SD, .SDcols=c(1,seq(9, 178, by = 4))] # get ncol(postZ), Prob(Z==2)
setkey(postZ_cons,V1)
colnames(postZ_cons)[2:44] <- species_name[1:43,1]

#### 1. get GC content for these CNEE compared with overall ####
# overall GC
gc_ratio0 = sum(GC[,-1])/sum(Tot[,-1]) # 0.3769491
gc_ratio0 = rowSums(GC[,-1])/rowSums(Tot[,-1])

gc_ratio_ratites0 = rowSums(GC[,ratites, with=F])/rowSums(Tot[,ratites, with=F])
gc_ratio_neo0 = rowSums(GC[,neognathae, with=F])/rowSums(Tot[,neognathae, with=F])

sel =  anno[,1]
gc_ratio1 = sum(GC[Elements %in% sel,-1])/sum(Tot[Elements %in% sel,-1]) #0.3521309

gc_ratio1 = rowSums(GC[Elements %in% sel,-1])/rowSums(Tot[Elements %in% sel,-1]) #
gc_ratio_ratites1 = rowSums(GC[Elements %in% sel,ratites, with=F])/rowSums(Tot[Elements %in% sel,ratites, with=F]) #
gc_ratio_neo1 = rowSums(GC[Elements %in% sel,neognathae, with=F])/rowSums(Tot[Elements %in% sel,neognathae, with=F]) #

x1 = c("All CNEEs", "HAR CNEEs")
x2 = c("all species",  "neognathaes","ratites")
label = paste(rep(x1, each=3), rep(x2,2), sep="\n")
pdf("Figures/GC/GC_HAR_boxplot.pdf", width=12)
boxplot(list(gc_ratio0, gc_ratio_neo0, gc_ratio_ratites0,gc_ratio1,gc_ratio_neo1, gc_ratio_ratites1), 
        names= label , outline = F, col = rep(3:4, each=3),notch=T, ylab = "GC content (%GC)" )
abline(h=0.3769491, col=2, lty=2 )
dev.off()

# tests
t.test(gc_ratio0, gc_ratio1, alternative = "less")
ks.test(gc_ratio0, gc_ratio1, alternative = "greater")


#### 2. Compare GC content for accelerated vs. conserved branches ####
GC_sel = rowSums(GC[sel,-1])/rowSums(Tot[sel,-1])
GC_sel1 = rowSums(GC[sel1,-1])/rowSums(Tot[sel1,-1])

# reorder columns of GC to be the same as postZ
GC1 = GC[, c("Elements", species_name[1:43,1]),with=F]
Tot1 = Tot[, c("Elements",species_name[1:43,1]),with=F]
setkey(GC1,Elements)
setkey(Tot1,Elements)

# get accelerated elements
sel <- score[which(score$loglik_Full - score$loglik_Null > 10),"ID"]
# get accelerated elements within HAR
sel1 <- intersect(anno[,1], sel) #100

#get accelerated species
selPostZ = (postZ_acc[sel,2:44, with=F] > 0.7)
GC_sel_acc = rowSums(GC1[sel,-1] * selPostZ)/rowSums(Tot1[sel,-1] * selPostZ)

# get conserved species
selPostZ = (postZ_cons[sel,2:44, with=F] > 0.7)
GC_sel_cons = rowSums(GC1[sel,-1] * selPostZ)/rowSums(Tot1[sel,-1] * selPostZ)

#get accelerated species for HAR
selPostZ = (postZ_acc[sel1,2:44, with=F] > 0.7)
GC_sel1_acc = rowSums(GC1[sel1,-1] * selPostZ)/rowSums(Tot1[sel1,-1] * selPostZ)

# get conserved species for HAR
selPostZ = (postZ_cons[sel1,2:44, with=F] > 0.7)
GC_sel1_cons = rowSums(GC1[sel1,-1] * selPostZ)/rowSums(Tot1[sel1,-1] * selPostZ)


x1 = c("32157 CNEEs", "100 HAR CNEEs")
x2 = c("all species",  "conserved species","accelerated species")
label = paste(rep(x1, each=3), rep(x2,2), sep="\n")

pdf("Figures/GC/GC_HAR_boxplot_acceleration.pdf", width=12)
boxplot(list(GC_sel, GC_sel_cons, GC_sel_acc,GC_sel1, GC_sel1_cons, GC_sel1_acc),
        names= label , outline = F, col = rep(3:4, each=3),notch=T, ylab = "GC content (%GC)", main="accelerated CNEEs")
abline(h=0.3769491, col=2, lty=2 )
abline(h=0.3521309, col=5, lty=2 )
legend("topleft", c("all CNEEs", "HAR CNEEs"),title="mean GC", lty=2, col=c(2,5))
dev.off()

# test
diff = GC_sel2_acc - GC_sel2_cons
t.test(diff[!is.na(diff)], alternative = "g")

#### 3 get ratite-accelerated elements ####
# sel2 <- intersect(anno[,1], score[which(score$logBF1>10 & score$logBF2 > 2),"ID"])
# 
# # get GC for the ratite-accelerated elements for accelerated lineages
# GC_sel2 = rowSums(GC[sel2,-1])/rowSums(Tot[sel2,-1])
# 
# # get accelerated species
# selPostZ = (postZ_acc[sel2,2:44, with=F] > 0.7)
# GC_sel2_acc = rowSums(GC1[sel2,-1] * selPostZ)/rowSums(Tot1[sel2,-1] * selPostZ)
# 
# # get conserved species
# selPostZ = (postZ_cons[sel2,2:44, with=F] > 0.7)
# GC_sel2_cons = rowSums(GC1[sel2,-1] * selPostZ)/rowSums(Tot1[sel2,-1] * selPostZ)
# 
# x1 = c("22436\naccelerated CNEEs", "820\nratite-acc CNEEs")
# x2 = c("all species",  "conserved species","accelerated species")
# label = paste(rep(x1, each=3), rep(x2,2), sep="\n")
# 
# pdf("Figures/GC/GC_boxplot_acceleration2.pdf", width=12)
# boxplot(list(GC_sel1, GC_sel1_cons, GC_sel1_acc,GC_sel2, GC_sel2_cons, GC_sel2_acc),
#         names= label , outline = F, col = rep(3:4, each=3),notch=F, ylab = "GC content (%GC)")
# abline(h=0.3521309, col=2, lty=2 )
# dev.off()
