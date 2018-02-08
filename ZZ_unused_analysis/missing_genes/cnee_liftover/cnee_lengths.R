setwd("~/Desktop/cnee/edited/")
path=("~/Desktop/cnee/edited/")

#file.names <- dir(path, pattern =".txt")

#read in all files that end in .txt from maker/out and place them into one dataframe
# cnee.out<-""
# for(i in 1:length(file.names)){
#   file <- read.delim(file.names[i])
#   cnee.out <- Reduce(function(x, y) merge(x, y, all=TRUE), list(cnee.out,file))
# }
# 
# drops <- c("x") #create a list of the columns that will be uninformative soon and we'd like to drop (just the extra x)
# cnee.out <- cnee.out[,!(names(cnee.out) %in% drops)]
# 
# write.csv(cnee.out, "liftover_cnee_lengths.csv", row.names = FALSE)
cnee.out <- read.csv("liftover_cnee_lengths.csv")
#this is too big, so subsetting.
#cnee.out$random <- sample(1:6,nrow(cnee.out))

# cnee.out1 <- subset(cnee.out, cnee.out$random==1)
# cnee.out2 <- subset(cnee.out, cnee.out$random==2)
# cnee.out3 <- subset(cnee.out, cnee.out$random==3)
# cnee.out4 <- subset(cnee.out, cnee.out$random==4)
# cnee.out5 <- subset(cnee.out, cnee.out$random==5)
# cnee.out6 <- subset(cnee.out, cnee.out$random==6)
# 
# nrow(cnee.out1)+nrow(cnee.out2)+nrow(cnee.out3)+nrow(cnee.out4)+nrow(cnee.out5)+nrow(cnee.out6)
# write.csv(cnee.out1, "liftover_cnee_lengths1.csv", row.names = FALSE)
# write.csv(cnee.out2, "liftover_cnee_lengths2.csv", row.names = FALSE)
# write.csv(cnee.out3, "liftover_cnee_lengths3.csv", row.names = FALSE)
# write.csv(cnee.out4, "liftover_cnee_lengths4.csv", row.names = FALSE)
# write.csv(cnee.out5, "liftover_cnee_lengths5.csv", row.names = FALSE)
# write.csv(cnee.out6, "liftover_cnee_lengths6.csv", row.names = FALSE)


#sanity check - nothing should be longer than galGal - we hope.
#allanimals <- c('allMis','anoCar','casCas','galGal')
allanimals<-c('allMis','allSin','anaPla','anoCar','aptFor','aptHaa','aptOwe','aptRow','balReg','calAnn','casCas','chaPel','chaVoc','cheMyd','chrPic','colLiv','corBra','croPor','cryCin','cucCan','droNov','eudEle','falPer','ficAlb','fulGla','galGal','gavGan','halLeu','lepDis','melGal','melUnd','mesUni','nipNip','notPer','picPub','pseHum','pygAde','rheAme','rhePen','strCam','taeGut','tinGut')
cnee.out$max_test <- apply(cnee.out[,allanimals], 1, function(x) max(x, na.rm=T))
all(cnee.out$galGal==cnee.out$max_test) #TRUE

colnames(cnee.out)
ratite <- c("aptHaa","aptOwe","aptRow","casCas","droNov","rheAme","rhePen","strCam")
flying <-c("galGal","anaPla","aptFor","balReg","calAnn","chaPel","chaVoc","colLiv","corBra","cryCin","cucCan","eudEle","falPer","ficAlb","fulGla","halLeu","lepDis","melGal","melUnd","mesUni","nipNip","notPer","picPub","pseHum","pygAde","taeGut","tinGut")
herps <- c("allMis","allSin","cheMyd","anoCar","chrPic","gavGan","croPor")
tinamou <- c("cryCin","eudEle","notPer","tinGut")

cnee.out$r_low <- apply(cnee.out[,ratite], 1, function(x) min(x, na.rm=T))
cnee.out$r_high <- apply(cnee.out[,ratite], 1, function(x) max(x, na.rm=T))
cnee.out$f_low <- apply(cnee.out[,flying], 1, function(x) min(x, na.rm=T))
cnee.out$t_low <- apply(cnee.out[,tinamou], 1, function(x) min(x, na.rm=T))
cnee.out$h_high <- apply(cnee.out[,herps], 1, function(x) max(x, na.rm=T))
cnee.out$h_low <- apply(cnee.out[,herps], 1, function(x) min(x, na.rm=F))

cnee.out$r_present <- (cnee.out$r_high/cnee.out$galGal)
cnee.out$f_present <- (cnee.out$f_low/cnee.out$galGal)
cnee.out$t_present <- (cnee.out$t_low/cnee.out$galGal)

ratite.lowfruit <- subset(cnee.out, cnee.out$r_present < 0.5 & cnee.out$r_present > -Inf & cnee.out$t_present > 0.95 & cnee.out$t_present < Inf & cnee.out$f_present > 0.95 & cnee.out$f_present < Inf)

na.lowfruit <- subset(cnee.out, cnee.out$h_high == -Inf & cnee.out$r_present < 0.5 & cnee.out$r_present > -Inf & cnee.out$t_present > 0.95 & cnee.out$t_present < Inf & cnee.out$f_present > 0.95 & cnee.out$f_present < Inf)

accel <- read.table("ratite_accel.tsv")

or.cnee.out <- cnee.out[with(cnee.out, order(ce)), ]
or.accel <- accel[with(accel, order(V1)), ]

#sanity checks
all(accel$V1==cnee.out$ce)
all(sort(accel$V1)==sort(cnee.out$ce))
all(or.accel$V1==or.cnee.out$ce)

#make accel column
or.cnee.out$accel <- or.accel$V2
#make fdr column
or.cnee.out$fdr_accel <- p.adjust(or.cnee.out$accel, method = "fdr",n=length(or.cnee.out$accel))

or.ratite.lowfruit <- subset(or.cnee.out, or.cnee.out$r_present < 0.5 & or.cnee.out$r_present > -Inf & or.cnee.out$t_present > 0.95 & or.cnee.out$t_present < Inf & or.cnee.out$f_present > 0.95 & or.cnee.out$f_present < Inf)

or.ratite.lowfruit <- subset(or.cnee.out, or.cnee.out$r_present < 0.5 & or.cnee.out$r_present > -Inf & or.cnee.out$t_present > 0.95 & or.cnee.out$t_present < Inf & or.cnee.out$f_present > 0.95 & or.cnee.out$f_present < Inf)

or.ratite.lowfruit.accel <- subset(or.cnee.out, or.cnee.out$r_present < 0.5 & or.cnee.out$r_present > -Inf & or.cnee.out$t_present > 0.95 & or.cnee.out$t_present < Inf & or.cnee.out$f_present > 0.95 & or.cnee.out$f_present < Inf & or.cnee.out$fdr_accel < 0 & or.cnee.out$fdr_accel > -1e-02)


#11075 obs of convergent ratite acceleration in ce's greater than 100 bp in chicken
or.ratite.accel <- subset(or.cnee.out, or.cnee.out$galGal > 100 & or.cnee.out$fdr_accel < 0 & or.cnee.out$fdr_accel > -1e-02)
#2869 obs of convergent ratite deceleration in ce's greater than 100 bp in chicken
or.ratite.decel <- subset(or.cnee.out, or.cnee.out$galGal > 100 & or.cnee.out$fdr_accel > 0 & or.cnee.out$fdr_accel < 1e-02)
hist(or.ratite.decel$galGal,or.ratite.accel$galGal)

annot <- read.table("ce_annotation.tsv",header = T)
all(annot$id==or.cnee.out$ce) #TRUE!
#add in annotation column
or.cnee.out$annot <- annot$gene

library(ggplot2)
#this will crash everything!
#ggp <- ggplot(data.frame(working),aes(x=working$annot))
# # counts
#ggp + geom_histogram(fill="lightgreen")


install.packages("data.table")        # install it
library(data.table) #allows manipulation
dt <- data.table(or.cnee.out) #create data.table of our working dataframe
dt_sig <- subset(dt, fdr_accel < 0 & or.cnee.out$fdr_accel > -1e-02) #and create one that only includes significant values

logplot.data <- dt[,.N,by=annot] #count all CNEEs for each exon
logplot.sig.data <- dt_sig[,.N,by=annot] #count all sig CNEEs for each exon
#order the datasets
logplot.data <- logplot.data[with(logplot.data, order(logplot.data$annot)), ]
logplot.sig.data <- logplot.sig.data[with(logplot.sig.data, order(logplot.sig.data$annot)), ]

#rename the "N" column from sig data to "sig"
library(plyr)
logplot.sig.data <- rename(logplot.sig.data,c("N"="sig"))

#combine the datasets
logplot.master <- Reduce(function(x, y) merge(x, y, all=TRUE, by="annot"), list(logplot.data,logplot.sig.data))
#there are no 0's, so replace NA with 0
logplot.master[is.na(logplot.master)] <- 0
logplot.master$percent <- logplot.master$sig/logplot.master$N

#basic scatter
p <- ggplot(logplot.master, aes(sig, N))
p + geom_point(position='jitter')

#log the y
p <- ggplot(logplot.master, aes(sig, N))
p + geom_point(position='jitter') + coord_trans(y = "log10")

#playing with point labels
p <- ggplot(logplot.master, aes(x=sig, y= N, label=annot))
p + geom_point(position='jitter') + coord_trans(y = "log10") + geom_text(aes(label=ifelse(sig>150,as.character(annot),'')),hjust=1, vjust=0)

#the final CE plot (without splitting CEs with multiple associated genes)
p <- ggplot(logplot.master, aes(x=sig, y= N, label=annot))
p + geom_point(position='jitter',alpha=0.6,size=1, color=ifelse(logplot.master$percent>=0.2,"magenta","white")) + coord_trans(y = "log10") + theme(panel.background = element_rect(fill = "black"), plot.title = element_text(size=10)) + ylab("Log Number of CEs associated with gene") + xlab("Number of CEs accelerated in ratites for gene (fdr corrected p value > -0.01)") + labs(title = "Each gene's number of associated CEs versus the number accelerated in ratite lineages. Pink points have > 0.2 total CEs accelerated in ratites", size=2)

#going back to subset only the CNEEs before carrying out same analyses

or.cnee.out$class <- annot$class
unique(or.cnee.out$class) #to go from CEs (everything) to CNEEs, we need only genic_non_exonic and intergenic
true.CNEEs <- subset(or.cnee.out, or.cnee.out$class %in% c("genic_non_exonic","intergenic"))
unique(true.CNEEs$class)
ct <- data.table(true.CNEEs) #create data.table
ct_sig <- subset(ct, fdr_accel < 0 & fdr_accel > -1e-02) #create sig data.table

test <- ct[,.N,by=class] #sanity check - only have two classes.
rm(test)

# total CNEE obs is over 1.7 M.  
# 43210 obs of convergent ratite acceleration in CNEEs greater than 100 bp in chicken (2.5%)
accel.ct <- subset(ct, ct$fdr_accel < 0 & ct$fdr_accel > -1e-02)
# 4433 obs of convergent ratite deceleration in CNEEs greater than 100 bp in chicken (0.03%)
decel.ct <- subset(ct, ct$fdr_accel > 0 & ct$fdr_accel < 1e-02)

# total over 100 bp in chicken is 110274.
onehunplus.ct <- subset(ct, ct$galGal > 100)
# 9930 obs of convergent ratite acceleration in CNEEs greater than 100 bp in chicken (9%)
ohp.accel.ct <- subset(ct, ct$galGal > 100 & ct$fdr_accel < 0 & ct$fdr_accel > -1e-02)
# 1881 obs of convergent ratite deceleration in CNEEs greater than 100 bp in chicken (1%)
ohp.decel.ct <- subset(ct, ct$galGal > 100 & ct$fdr_accel > 0 & ct$fdr_accel < 1e-02)

logplot.ct.data <- ct[,.N,by=annot] #count all CNEEs for each exon
logplot.ct.sig.data <- ct_sig[,.N,by=annot] #count all sig CNEEs for each exon
#order the datasets
logplot.ct.data <- logplot.ct.data[with(logplot.ct.data, order(logplot.ct.data$annot)), ]
logplot.ct.sig.data <- logplot.ct.sig.data[with(logplot.ct.sig.data, order(logplot.ct.sig.data$annot)), ]

#combine the datasets
logplot.ct.master <- Reduce(function(x, y) merge(x, y, all=TRUE, by="annot"), list(logplot.ct.data,logplot.ct.sig.data))
#there are no 0's, so replace NA with 0
logplot.ct.master[is.na(logplot.ct.master)] <- 0
#add the percent column
logplot.ct.master$percent <- logplot.ct.master$N.y/logplot.ct.master$N.x

#the final CNEE plot (without splitting CNEEs with multiple associated genes)
p.ct <- ggplot(logplot.ct.master, aes(x=N.y, y= N.x, label=annot))
p.ct + geom_point(position='jitter',alpha=0.6,size=1, color=ifelse(logplot.ct.master$percent>=0.2,"magenta","white")) + coord_trans(y = "log10") + theme(panel.background = element_rect(fill = "black"), plot.title = element_text(size=10)) + ylab("Log Number of CNEEs associated with gene") + xlab("Number of CNEEs accelerated in ratites for gene (fdr corrected p value > -0.01)") + labs(title = "Each gene's number of associated CNEEs versus the number accelerated in ratite lineages. Pink points have > 0.2 total CNEEs accelerated in ratites", size=2)


NKX25 <- subset(true.CNEEs, true.CNEEs$annot =="CGNC:2083,GeneID:396073")
sigNXK <- subset(NKX25, fdr_accel < 0 & fdr_accel > -1e-02)
loNKX <- subset(NKX25, NKX25$r_present < 0.5 & NKX25$r_present > -Inf & NKX25$t_present > 0.95 & NKX25$t_present < Inf & NKX25$f_present > 0.95 & NKX25$f_present < Inf)
loNKX <- subset(NKX25, NKX25$r_low < 0.5*NKX25$galGal & NKX25$r_present > -Inf & NKX25$t_present > 0.95 & NKX25$t_present < Inf & NKX25$f_present > 0.95 & NKX25$f_present < Inf)


SHOX2 <- true.CNEEs[grep("GeneID:777244", true.CNEEs$annot), ]
sigSHOX2 <- subset(SHOX2, fdr_accel < 0 & fdr_accel > -1e-02)
loSHOX2 <- subset(SHOX2, SHOX2$r_present < 0.5 & SHOX2$r_present > -Inf & SHOX2$t_present > 0.95 & SHOX2$t_present < Inf & SHOX2$f_present > 0.95 & SHOX2$f_present < Inf)
loSHOX2 <- subset(SHOX2, SHOX2$r_low < 0.5*SHOX2$galGal & SHOX2$r_present > -Inf & SHOX2$t_present > 0.95 & SHOX2$t_present < Inf & SHOX2$f_present > 0.95 & SHOX2$f_present < Inf)

unbiased.both.sets <- subset(true.CNEEs, true.CNEEs$r_present < 0.5 & true.CNEEs$r_present > -Inf & true.CNEEs$t_present > 0.95 & true.CNEEs$t_present < Inf & true.CNEEs$f_present > 0.95 & true.CNEEs$f_present < Inf & true.CNEEs$fdr_accel < 0.01 & true.CNEEs$accel < 0)

true.CNEEs$abs_fdr <- p.adjust(abs(true.CNEEs$accel), method = "fdr",n=length(true.CNEEs$accel))

write.csv(true.CNEEs, "phil_cnee_table.csv", row.names = FALSE)
sample <- true.CNEEs[sample(rownames(true.CNEEs),n=5000),]
plot(sample$abs_fdr,sample$fdr_accel)
# ggplot()+aes(or.cnee.out$galGal)+geom_histogram(colour="lightblue", fill="black", binwidth=25)
# ggplot()+aes(or.cnee.out$galGal)+geom_histogram(colour="lightblue", fill="black", binwidth=1)+ xlim(0,50)
# ggplot()+aes(or.cnee.out$galGal)+geom_histogram(colour="lightblue", fill="black", binwidth=125)+ xlim(1000,3000)
# ggplot()+aes(or.cnee.out$galGal)+geom_histogram(colour="lightblue", fill="black", binwidth=3)+ xlim(200,1000)
# 
# ggplot()+aes(or.ratite.accel$galGal)+geom_histogram(colour="pink", fill="white")
# ggplot()+aes(or.ratite.decel$galGal)+geom_histogram(colour="blue", fill="white")
# drops <- c("r_perc_loss") #create a list of the columns that will be uninformative soon and we'd like to drop (just the extra x)
# cnee.out <- cnee.out[,!(names(cnee.out) %in% drops)]
# 
# hist(cnee.out$r_perc_loss)

or.cnee.out <- read.csv("phil_cnee_table.csv")


