setwd("~/Desktop/orf_missing/")
listy <- list.files(path = ".")

rA <- read.delim("rheAme_orfs.txt")
cP <- read.delim("croPor_orfs.txt")
nP <- read.delim("notPer_orfs.txt")
fP <- read.delim("falPer_orfs.txt")
mU <- read.delim("mesUni_orfs.txt")

all.equal(rA$ChickenSeqID,cP$ChickenSeqID)
all.equal(rA$ChickenSeqID,nP$ChickenSeqID)
all.equal(rA$ChickenSeqID,fP$ChickenSeqID)
all.equal(rA$ChickenSeqID,mU$ChickenSeqID)

path = "~/Desktop/orf_missing/"
out.file<-""
file.names <- dir(path, pattern =".txt")
for(i in 1:length(file.names)){
  file <- read.delim(file.names[i])
  out.file <- Reduce(function(x, y) merge(x, y, all=TRUE), list(out.file,file))
}
write.table(out.file, file = "testing.txt",sep="\t")

out.file[grep("100170842,Genbank:NP_001124216.1",out.file$ChickenSeqID),]
rA[grep("100170842,Genbank:NP_001124216.1",rA$ChickenSeqID),]



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


