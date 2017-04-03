#Loading the bayesian CNEE files and writing out bed files with genomic locations

setwd("~/Dropbox/Doctorate/atac/ChickATAC/")
by.just.conAcc<- read.table("final_bayesian_conv_Jan2017.txt",header=T) 
cnee <- read.table("final_cnees_long.bed")
final_bayes <- merge(cnee, by.just.conAcc, by.x='V4', by.y='Row.names')

table(final_bayes$clades) #check that these two are the same
table(by.just.conAcc$clades)

looseBioGeo <- (final_bayes[which(final_bayes$clades>=2),])
strictBioGeo <- (final_bayes[which(final_bayes$clades>=4),])
looseDollo <- (final_bayes[which(final_bayes$clades>=2 & final_bayes$anoDid == 1 | final_bayes$strCam == 1),])
strictDollo <- (final_bayes[which(final_bayes$clades>=2 & final_bayes$anoDid == 1 & final_bayes$strCam == 1),])

write.table(looseBioGeo,sep = "\t",file = "looseBioGeo.bed",quote = FALSE,col.names = F,row.names = F)
write.table(strictBioGeo,sep = "\t",file = "strictBioGeo.bed",quote = FALSE,col.names = F,row.names = F)      
write.table(looseDollo,sep = "\t",file = "looseDollo.bed",quote = FALSE,col.names = F,row.names = F)
write.table(strictDollo,sep = "\t",file = "strictDollo.bed",quote = FALSE,col.names = F,row.names = F)



#Following bedtools overlap, and python permutations, load the outfiles in

setwd("~/Dropbox/Doctorate/atac/MarchPeakAnalysis/NewBayesDistributions/")

listDir <- list.files(path=getwd(), pattern=".txt", full.names=FALSE)
DistList <- lapply(listDir, function(x) read.table(x, sep="\t",header=F)) #read the files in into a list
for (i in 1:length(listDir)){
  listDir[[i]] <- paste(strsplit(listDir[[i]], "[.]")[[1]][1:3],collapse="")
} #simplify names
names(DistList)<-listDir #associate names

for (i in 1:length(DistList)){
  x <- (eval(parse(text=(paste0('DistList$',listDir[[i]])))))
  x$Library <- listDir[[i]]
  DistList[[i]] <- x
} #make a new column called library that grabs the dataframe name

#make corrected column (overlapping bases/total bases in permutation)
DistList <- lapply(DistList, function(x) {
  x$corrected <- x$V2/x$V1
  return(x)
})


#bind the appropriate libraries (1+2,3+4,5+6, etc.)
for (i in 1:length(DistList)) {
  if (!i %% 2){
    next
  }
  assign(paste0("DF", i), do.call("rbind", DistList[i:(i+1)]))
}

#need a new list of DFs
DFs <- list(DF1,DF3,DF5,DF7,DF9,DF11,DF13,DF15,DF17,DF19,DF21,DF23,DF25,DF27,DF29,DF31,DF33,DF35,DF37,DF39,DF41,DF43,DF45,DF47,DF49,DF51,DF53,DF55,DF57,DF59,DF61,DF63)

#clean up the unlisted dataframes
rm(DF1,DF3,DF5,DF7,DF9,DF11,DF13,DF15,DF17,DF19,DF21,DF23,DF25,DF27,DF29,DF31,DF33,DF35,DF37,DF39,DF41,DF43,DF45,DF47,DF49,DF51,DF53,DF55,DF57,DF59,DF61,DF63)

#we have merged the correct libraries, so now we want to rename the dataframes correctly
newList <- list()
for (i in 1:length(listDir)) {
  if (!i %% 2){
    next
  }
  newList[[i]] <- paste(strsplit(listDir[[i]], "[_]")[[1]][4:5],collapse="")
}

cleanNewList<-newList[!sapply(newList, is.null)] #remove nulls from list
names(DFs) <- cleanNewList #associate the right names!


### Loading in the true wo overlap for the libraries ###
#note that this requires the python script March_Peak_trueIntersectWo_editor.py to go from the commented out file below to the ..._DFnames.txt that will work with this code.

setwd("~/Dropbox/Doctorate/atac/MarchPeakAnalysis/")
#newBayes.observed <- read.table(file = "woForNewBayesCategories.txt", sep=" ",header=F)
newBayes.observed <- read.table(file = "woForNewBayesCategories_DFnames.txt", sep="\t",header=F)
newBayes.observed$corrected  <- newBayes.observed$V3/newBayes.observed$V2
row.names(newBayes.observed) <- newBayes.observed$V1 #assign column values to row names

addObs <- function(i) {
  DFs[[i]]$trueBayes <- newBayes.observed[names(DFs)[i],]$corrected
  return(DFs[[i]])
} #define a function that adds a "trueBayes" column to each DF in DFs, generating the value from the $corrected column in the row in newBayes.observed that has the same name as the dataframe 

DFs <- lapply(seq_along(DFs),addObs) #apply the function
names(DFs) <- cleanNewList #name the dataframes again

### PLOTTING ###

#this works to plot here, but we also want to write out
# lapply(DFs,function(x){
#   ggplot(x, aes(x$corrected, fill = x$Library)) + geom_density(alpha = 0.2) + theme(legend.position="bottom")
# })

#now we want to create a plot for each DF in DFs and name it appropriately
setwd("~/Dropbox/Doctorate/atac/MarchPeakAnalysis/NewBayesDistributions/pdfOut/")

for(y in 1:length(DFs)) {
  DFname<-names(DFs[y])
  pdf(paste("plot", DFname, ".pdf", sep = ""))
  DFvariable<-(paste0("DFs$",DFname))
  DFtrueBayes<-(paste0("DFs$",DFname,"$trueBayes"))
  x<-ggplot(eval(parse(text = DFvariable)), aes(corrected, fill = Library)) + geom_density(alpha = 0.2) + theme(legend.position="bottom") + geom_vline(xintercept=eval(parse(text = DFtrueBayes)))
  print(x)
  dev.off()
}

### Permutation test p-values ###

#Make a column that is called GeneralLib denoting whether the row falls into the distribution of CNEEs sampled from "All" or "Accel" based on the number of CNEEs in the initial Bayes bed file
DFs <- lapply(DFs, function(x) {
  x$GeneralLib <- "All"
  x[100001:200000,"GeneralLib"] <- "Accel" #change to Accel for bayes sampled from bayes
  x$mean <- mean(x$corrected[1:100000]) #take the mean of the All 
  x[100001:200000,"mean"]<-mean(x$corrected[100001:200000]) #change the bayes to bayes mean
  x$difMean <- x[100001,"mean"]-x[1,"mean"] #take the diff of the mean
  x$countHigherDist <- sum(x$corrected[1:100000] >= x[100001,"mean"])
  x$countHigherTrue <- sum(x$corrected[1:100000] >= x$trueBayes)
  x$p.HigherDist <- x$countHigherDist/100000
  x$p.HigherTrue <- x$countHigherTrue/100000
  return(x)
})

#test that this worked correctly (it did):
# lapply(DFs, function(x) {
#   table(x$Library,x$GeneralLib )
# })

lapply(DFs, function(x) {
  table(x$Library,x$p.HigherTrue)
})

lapply(DFs, function(x) {
  table(x$Library,x$p.HigherDist)
})

lapply(DFs, function(x) {
  table(x$countHigherDist,x$p.HigherDist)
})

lapply(DFs, function(x){
  #print(names(x))
  #x[1,"p.HigherDist"]
  x[1,"p.HigherTrue"]
})