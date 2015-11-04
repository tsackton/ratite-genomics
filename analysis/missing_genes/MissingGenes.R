setwd("~/Desktop/out_of_stops/")
path = "~/Desktop/out_of_stops/"

rA <- read.delim("rheAme_stops.txt")
rA$best <- apply(rA[,c(2,3,4)],1,min)
rA$Gene <- sub(",.*","",rA$ChickenSeqID,perl=TRUE)
rasub <- subset(rA, select=c(6,7))
rasub <- rasub[order(rasub$Gene,rasub$best),]
ragene <- rasub[!duplicated(rasub$Gene),]
table(ragene$best)

best.out<-""
file.names <- dir(path, pattern =".txt")
for(i in 1:length(file.names)){
  file <- read.delim(file.names[i])
  file$best <- apply(file[,c(2,3,4)],1,min)
  file$Gene <- sub(",.*","",file$ChickenSeqID,perl=TRUE)
  filesub <- subset(file, select=c(6,7))
  filesub <- filesub[order(filesub$Gene,filesub$best),]
  gene <- filesub[!duplicated(filesub$Gene),]
  best.out <- Reduce(function(x, y) merge(x, y, all=TRUE), list(best.out,gene))
}



out.file<-""
file.names <- dir(path, pattern =".txt")
for(i in 1:length(file.names)){
  file <- read.delim(file.names[i])
  out.file <- Reduce(function(x, y) merge(x, y, all=TRUE), list(out.file,file))
}

drops <- c("x") #create a list of the columns that will be uninformative soon and we'd like to drop
orf.table <- out.file[,!(names(out.file) %in% drops)]

write.table(out.file, file = "testing.txt",sep="\t")
#rA <- read.delim("rheAme_orfs.txt")
#cP <- read.delim("croPor_orfs.txt")
#nP <- read.delim("notPer_orfs.txt")
#fP <- read.delim("falPer_orfs.txt")
#mU <- read.delim("mesUni_orfs.txt")

#all.equal(rA$ChickenSeqID,cP$ChickenSeqID)
#all.equal(rA$ChickenSeqID,nP$ChickenSeqID)
#all.equal(rA$ChickenSeqID,fP$ChickenSeqID)
#all.equal(rA$ChickenSeqID,mU$ChickenSeqID)

#ratite.lowfruit <- writer_names_df <- subset(orf.table, orf.table$rheAme > 2 & orf.table$rhePen > 2 & orf.table$casCas > 2 & orf.table$droNov > 2 & orf.table$strCam > 2 & orf.table$aptFor > 2 & orf.table$aptHaa > 2 & orf.table$aptOwe > 2 & orf.table$aptRow > 2 & orf.table$eudEle <= 2 & orf.table$cryCin <= 2 & orf.table$notPer <= 2 & orf.table$tinGut <= 2)

ratite.lowfruit <- writer_names_df <- subset(orf.table, orf.table$rheAme > 2 & orf.table$rhePen > 2 & orf.table$casCas > 2 & orf.table$droNov > 2 & orf.table$strCam > 2 & orf.table$aptHaa > 2 & orf.table$aptOwe > 2 & orf.table$aptRow > 2 & orf.table$eudEle <= 2 & orf.table$cryCin <= 2 & orf.table$notPer <= 2 & orf.table$melUnd <= 2 & orf.table$pseHum <= 2 & orf.table$colLiv <= 2 & orf.table$falPer <= 2 & orf.table$anaPla <= 2 & orf.table$fulGla <= 2 & orf.table$lepDis <= 2 & orf.table$corBra <= 2 & orf.table$mesUni <= 2 & orf.table$picPub <= 2  & orf.table$calAnn <= 2 & orf.table$pygAde <= 2 & orf.table$aptFor <= 2 & orf.table$chaVoc <= 2 & orf.table$nipNip <= 2 & orf.table$cucCan <= 2 & orf.table$balReg <= 2 & orf.table$halLeu <= 2 & orf.table$chaPel & orf.table$tinGut <= 2 )

#ratite_indexes = c(1,2,3,4)
#fly_indexes= same but diff
orf.table$badRats <- apply(orf.table[ratite_indexes,], 1, function(x) sum(x > 2)/length(x))

#compute new column apply(orf.table[ratite_indexes,], 1, function(x) sum(x > 2)/length(x))

new.lowfruit 

out.file[grep("100170842,Genbank:NP_001124216.1",out.file$ChickenSeqID),]
rA[grep("100170842,Genbank:NP_001124216.1",rA$ChickenSeqID),]
rAsort <- rA[order(rA$ChickenSeqID),]


Reduce(function(x, y) merge(x, y, all=TRUE), list(rA,cP,nP,fP,mU))

listy <-list.files(path=".")
LL <- list()
for(file in listy)
  LL <- c(LL,c=strsplit(file,"_")[[1]][1])

for(name in LL)
  this <- name
  this <- read.delim(paste(name, "_orfs.txt", sep=""))

strsplit(list.files(path=".")[1],"_")[[1]][1]

ldf <- list() # creates a list
rA <- read.delim("rheAme_orfs.txt")
listcsv <- dir(pattern = "*.txt") # creates the list of all the csv files in the directory
for (k in 1:length(listcsv)){
  ldf[[k]] <- read.delim(listcsv[k])
}
str(ldf[[1]]) 



#All_things <- rA
#All_things$croPor <- cP$croPor
#listy <- list.files(path = ".")

