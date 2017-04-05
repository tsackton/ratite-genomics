#Loading the bayesian CNEE files and writing out bed files with genomic locations

setwd("~/Dropbox/Doctorate/atac/ChickATAC/")
by.just.conAcc<- read.table("final_bayesian_conv_Jan2017.txt",header=T) 
cnee <- read.table("final_cnees_long.bed")
final_bayes <- merge(cnee, by.just.conAcc, by.x='V4', by.y='Row.names')

table(final_bayes$clades) #check that these two are the same
table(by.just.conAcc$clades)

looseBioGeo <- (final_bayes[which(final_bayes$clades>=2),])
strictBioGeo <- (final_bayes[which(final_bayes$clades>=4),])
looseDollo <- (final_bayes[which(final_bayes$clades>=2 & (final_bayes$anoDid == 1 | final_bayes$strCam == 1)),])
strictDollo <- (final_bayes[which(final_bayes$clades>=2 & final_bayes$anoDid == 1 & final_bayes$strCam == 1),])

write.table(looseBioGeo,sep = "\t",file = "looseBioGeo.bed",quote = FALSE,col.names = F,row.names = F)
write.table(strictBioGeo,sep = "\t",file = "strictBioGeo.bed",quote = FALSE,col.names = F,row.names = F)      
write.table(looseDollo,sep = "\t",file = "looseDollo.bed",quote = FALSE,col.names = F,row.names = F)
write.table(strictDollo,sep = "\t",file = "strictDollo.bed",quote = FALSE,col.names = F,row.names = F)



#Following bedtools overlap, and python permutations, load the outfiles in
#I have fixed the names of these files now, so there I have commented out the old way of dealing with them

#setwd("~/Dropbox/Doctorate/atac/MarchPeakAnalysis/NewBayesDistributions/")
setwd("~/Dropbox/Doctorate/atac/MarchPeakAnalysis/NewBayesDistributions/fixed_names/")

listDir <- list.files(path=getwd(), pattern=".txt", full.names=FALSE)
DistList <- lapply(listDir, function(x) read.table(x, sep="\t",header=F)) #read the files in into a list
# for (i in 1:length(listDir)){
#   listDir[[i]] <- paste(strsplit(listDir[[i]], "[.]")[[1]][1:3],collapse="")
# } #simplify names

for (i in 1:length(listDir)){
   listDir[[i]] <- paste(strsplit(listDir[[i]], "[.]")[[1]][1],collapse="")
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
# newList <- list()
# for (i in 1:length(listDir)) {
#   if (!i %% 2){
#     next
#   }
#   newList[[i]] <- paste(strsplit(listDir[[i]], "[_]")[[1]][4:5],collapse="")
# }

newList <- list()
for (i in 1:length(listDir)) {
  if (!i %% 2){
    next
  }
  newList[[i]] <- paste(strsplit(listDir[[i]], "[_]")[[1]][1:2],collapse="_")
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

View(DFs$Keel10_looseBioGeo)
### PLOTTING ###

#this works to plot here, but we also want to write out
# lapply(DFs,function(x){
#   ggplot(x, aes(x$corrected, fill = x$Library)) + geom_density(alpha = 0.2) + theme(legend.position="bottom")
# })

#now we want to create a plot for each DF in DFs and name it appropriately
library(ggplot2)
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


DFs <- lapply(DFs, function(x) {
  x$GeneralLib <- "All" #Make a column that is called GeneralLib denoting whether the row falls into the distribution of CNEEs sampled from "All" or "Accel" based on the number of CNEEs in the initial Bayes bed file
  x[100001:200000,"GeneralLib"] <- "Accel" #change to Accel for bayes sampled from bayes
  x$mean <- mean(x$corrected[1:100000]) #take the mean of the All 
  x[100001:200000,"mean"]<-mean(x$corrected[100001:200000]) #change the bayes rows to bayes mean
  x$difMean <- x[100001,"mean"]-x[1,"mean"] #take the diff of the mean
  x$countHigherDist <- sum(x$corrected[1:100000] >= x[100001,"mean"]) #count the number of rows from the All condition that have values higher than the mean of the Bayes distribution
  x$countHigherTrue <- sum(x$corrected[1:100000] >= x$trueBayes[1:100000]) #count the number of rows from the All condition that have values higher than the true Bayes overlap
  x$p.HigherDist <- x$countHigherDist/100000 #divide the count by 100000
  x$p.HigherTrue <- x$countHigherTrue/100000 #divide the other count by 100000
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

### new comparative enrichment test ###

####SUMMARYSE from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper functions
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


#bind the appropriate libraries (1+2,3+4,5+6, etc.)
for (i in 1:length(DistList)) {
  if (!i %% 2){
    next
  }
  assign(paste0("DF", i), do.call("cbind", DistList[i:(i+1)]))
}

#need a new list of DFs
wideDFs <- list(DF1,DF3,DF5,DF7,DF9,DF11,DF13,DF15,DF17,DF19,DF21,DF23,DF25,DF27,DF29,DF31,DF33,DF35,DF37,DF39,DF41,DF43,DF45,DF47,DF49,DF51,DF53,DF55,DF57,DF59,DF61,DF63)
#clean up the unlisted dataframes
rm(DF1,DF3,DF5,DF7,DF9,DF11,DF13,DF15,DF17,DF19,DF21,DF23,DF25,DF27,DF29,DF31,DF33,DF35,DF37,DF39,DF41,DF43,DF45,DF47,DF49,DF51,DF53,DF55,DF57,DF59,DF61,DF63)

#made a function to change column names
ChangeCols <- function(x) {
  colnames(x)<-paste0("V",(1:ncol(x))) 
  return(x)
}

#and rename the columns
wideDFs <- lapply(wideDFs,ChangeCols)
names(wideDFs) <- cleanNewList

wideDFs <- lapply(wideDFs, function(x) {
  x$BayesMinusCNEE <- x$V8-x$V4
  x$basesPerKB <- (x$BayesMinusCNEE)*1000
  return(x)
})

#install.packages("tidyr")
#install.packages("dplyr")
library(tidyr)
library(dplyr)
library(ggplot2)

#with tidry and dplyr, we can do some nice new stuff.
#want to get gat values and associate them to main DF that we create below

setwd("~/Dropbox/Doctorate/atac/MarchPeakAnalysis/gatOut/")

listDir2 <- list.files(path=getwd(), pattern=".txt", full.names=FALSE)
DistList2 <- lapply(listDir2, function(x) read.table(x, sep="\t",header=T)) #read the files in into a list
for (i in 1:length(listDir2)){
  listDir2[[i]] <- paste(strsplit(listDir2[[i]], "[._]")[[1]][1],collapse="")
}  #simplify names
names(DistList2) <- listDir2

gatOut <- lapply(DistList2, function(x) {
  x <- x %>% tbl_df %>% select(annotation, observed, expected, fold, pvalue, qvalue)
  return(x)
}) #change each df in list to tbl and only keep columns of current interest


annotation_key = list("Bayes_S_D"="strictDollo", "Bayes_L_BG"="looseBioGeo", "Bayes_L_D"="looseDollo", "Bayes_S_BG"="strictBioGeo","BAYES"="oldBayes","PHAST"="oldPhast","CNEE"="AllCNEEs","TSS"="TSS") #necessary since the two dfs don't agree on what to call the different annotations

gatOut.wide <- do.call("rbind", gatOut)

widebind <- do.call("rbind", wideDFs)
sum.widebind <- summarySE(widebind, measurevar="basesPerKB",groupvars='V3')
wider.sum.widebind <- sum.widebind %>% tbl_df %>% separate(V3, c("Tissue","Bayes.List"), remove = T, extra = "drop")
BayesClasses <- wider.sum.widebind %>% distinct(Bayes.List)

merged_enrichment <- gatOut.wide %>% tibble::rownames_to_column() %>% separate(rowname, c("Tissue"), remove = T, extra = "drop") %>% mutate(bayes.key = unlist(annotation_key[as.character(annotation)])) %>% inner_join(wider.sum.widebind, ., by=c("Bayes.List" = "bayes.key", "Tissue" = "Tissue")) %>% mutate(pValue = ifelse(pvalue < 0.05, "p < 0.05", "p > 0.05"))


#par(mfrow=c(1,2)) this notation doesn't work with ggplot.  need to do the following
#install.packages("gridExtra")
library(gridExtra)

for (i in 1:length(BayesClasses$Bayes.List)){
  assign(paste0("p",i), merged_enrichment %>% tbl_df %>% filter(Bayes.List == BayesClasses$Bayes.List[[i]]) %>%
           ggplot(., aes(x=Tissue, y=basesPerKB, fill=pValue)) + geom_bar(stat="identity") +
           ggtitle(paste0("Differential Enrichment Across Tissue and Developmental Stage for ",BayesClasses$Bayes.List[[i]])) +
           ylab("Enrichment Per 1000 bases") + theme(plot.title = element_text(size = 9, face = "bold")) +
           scale_y_continuous(breaks = seq(-10,100, by=10), limits = c(-10,100)) +
           geom_errorbar(aes(ymin=(basesPerKB-ci), ymax=(basesPerKB+ci)),size=.3,width=.2,position=position_dodge(.9)) +
           scale_fill_manual(values = c("p > 0.05" = "darkgrey", "p < 0.05" = "springgreen4"))+
           theme(legend.position = "bottom"))
} 
grid.arrange(p1, p2, p3, p4, ncol=2) #plot these 4 plots in 2x2 figure.
