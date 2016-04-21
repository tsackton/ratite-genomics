setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/paml/status_check/")
status<-read.table("paml-status-20160412", header=F, stringsAsFactors=F)
names(status)<-c("hog", "model", "time")

status.rerun <- subset(status, time == "NOT_STARTED" | time == "RUNNING")
status.times <- subset(status, time != "NOT_STARTED" & time != "RUNNING")

boxplot(as.numeric(status.times$time)/60/60 ~ status.times$model, las=1, ylab="Hours Run Time", main="PAML Run Time")

with(droplevels(status[status$model!="br",]), barplot(table(model, time=="NOT_STARTED" | time == "RUNNING")[,1]/11271, ylim=c(0,1), las=1, ylab="Fraction Finished", main="Fraction Finished"))


with(droplevels(status[status$model!="br",]),table(model, time=="NOT_STARTED" | time == "RUNNING")[,2])

with(droplevels(status[status$model!="br",]),by(as.numeric(time)/60, model, mean, na.rm=T))*with(droplevels(status[status$model!="br",]),table(model, time=="NOT_STARTED" | time == "RUNNING")[,2])/60/24

#write out reruns
write.table(status.rerun[status.rerun$model=="ancrec","hog"], file="ancrec.rerun", sep="\t", row.names=F, col.names=F, quote=F)
