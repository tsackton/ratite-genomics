setwd("~/Desktop/cnee/edited/")
path=("~/Desktop/cnee/edited/")

file.names <- dir(path, pattern =".txt")

#read in all files that end in .txt from maker/out and place them into one dataframe
cnee.out<-""
for(i in 1:length(file.names)){
  file <- read.delim(file.names[i])
  cnee.out <- Reduce(function(x, y) merge(x, y, all=TRUE), list(cnee.out,file))
}

drops <- c("x") #create a list of the columns that will be uninformative soon and we'd like to drop (just the extra x)
cnee.out <- cnee.out[,!(names(cnee.out) %in% drops)]

write.csv(cnee.out, "liftover_cnee_lengths.csv", row.names = FALSE)
#this is too big, so subsetting.
cnee.out$random <- sample(1:6,nrow(cnee.out))

cnee.out1 <- subset(cnee.out, cnee.out$random==1)
cnee.out2 <- subset(cnee.out, cnee.out$random==2)
cnee.out3 <- subset(cnee.out, cnee.out$random==3)
cnee.out4 <- subset(cnee.out, cnee.out$random==4)
cnee.out5 <- subset(cnee.out, cnee.out$random==5)
cnee.out6 <- subset(cnee.out, cnee.out$random==6)

nrow(cnee.out1)+nrow(cnee.out2)+nrow(cnee.out3)+nrow(cnee.out4)+nrow(cnee.out5)+nrow(cnee.out6)
write.csv(cnee.out1, "liftover_cnee_lengths1.csv", row.names = FALSE)
write.csv(cnee.out2, "liftover_cnee_lengths2.csv", row.names = FALSE)
write.csv(cnee.out3, "liftover_cnee_lengths3.csv", row.names = FALSE)
write.csv(cnee.out4, "liftover_cnee_lengths4.csv", row.names = FALSE)
write.csv(cnee.out5, "liftover_cnee_lengths5.csv", row.names = FALSE)
write.csv(cnee.out6, "liftover_cnee_lengths6.csv", row.names = FALSE)


#sanity check - nothing should be longer than galGal - we hope.
#allanimals <- c('allMis','anoCar','casCas','galGal')
allanimals<-c('allMis','allSin','anaPla','anoCar','aptFor','aptHaa','aptOwe','aptRow','balReg','calAnn','casCas','chaPel','chaVoc','cheMyd','chrPic','colLiv','corBra','croPor','cryCin','cucCan','droNov','eudEle','falPer','ficAlb','fulGla','galGal','gavGan','halLeu','lepDis','melGal','melUnd','mesUni','nipNip','notPer','picPub','pseHum','pygAde','rheAme','rhePen','strCam','taeGut','tinGut')
cnee.out$max_test <- apply(cnee.out[,allanimals], 1, function(x) max(x, na.rm=T))
all(cnee.out$galGal==cnee.out$max_test)

allanimals<-c('allMis','allSin','anaPla','anoCar','aptFor','aptHaa','aptOwe','aptRow','balReg','calAnn','casCas','chaPel','chaVoc','cheMyd','chrPic','colLiv','corBra','croPor','cryCin','cucCan','droNov','eudEle','falPer','ficAlb','fulGla','galGal','gavGan','halLeu','lepDis','melGal','melUnd','mesUni','nipNip','notPer','picPub','pseHum','pygAde','rheAme','rhePen','strCam','taeGut','tinGut')