setwd("~/Desktop/")
rA <- read.delim("rheAme_orfs.txt",row.names = 1)
cP <- read.delim("croPor_orfs.txt",row.names = 1)
nP <- read.delim("notPer_orfs.txt",row.names = 1)
fP <- read.delim("falPer_orfs.txt",row.names = 1)
mU <- read.delim("mesUni_orfs.txt",row.names = 1)

all.equal(rA$ChickenSeqID,cP$ChickenSeqID)
all.equal(rA$ChickenSeqID,nP$ChickenSeqID)
all.equal(rA$ChickenSeqID,fP$ChickenSeqID)
all.equal(rA$ChickenSeqID,mU$ChickenSeqID)

All_things <- merge(rA,cP, by.x=row.names,by.y=row.names)
?merge
All_things <- rA
All_things$croPor <- cP$croPor
