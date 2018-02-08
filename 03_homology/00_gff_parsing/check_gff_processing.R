#read in all_info.tab

all<-read.table("gffs/all_info.summary", header=T, comment.char="", sep="\t", stringsAsFactors=F)
names(all)[1]="species"
all=subset(all, species != "#species")
table(all$species, all$biotype, useNA="ifany")

barplot(table(all$species, all$biotype, useNA="ifany")[,6], las=2)
