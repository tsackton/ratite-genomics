library(tidyr)
library(dplyr)
setwd("~/Projects/birds/genome_quality/final/R")
cont<-read.csv("../contiguity/allstats.final.csv", header=T, stringsAsFactors = F)

#old code not worth rewriting as dplyr
#clean up
cont$Assembly = sub(".fas.gz", "", cont$Assembly, fixed=T)
cont$Assembly = sub(".fa", "", cont$Assembly, fixed=T)
#ugly code to fix names
names(cont)=c("assembly", "num.scaf", "tot.size", "useful.size", "longest.scaf", "shortest.scaf", "num.scaf.gt.1k", "perc.scaf.gt.1k", "num.scaf.gt.10k", "perc.scaf.gt.10k", "num.scaf.gt.100k", "perc.scaf.gt.100k", "num.scaf.gt.1m", "perc.scaf.gt.1m", "num.scaf.gt.10m", "perc.scaf.gt.10m", "scaf.mean.size", "scaf.median.size", "scaf.n50", "scaf.l50", "scaf.percA", "scaf.percC", "scaf.percG", "scaf.percT", "scaf.percN", "scaf.perc.nonATGCN", "scaf.num.nonATGCN", "perc.scaffolded.ct", "perc.unscaffolded.ct", "avg.ct.per.scaf", "avg.N.bw.ct", "num.ct", "num.ct.scaffolded", "num.ct.unscaffolded", "ct.tot.size", "longest.ct", "shortest.ct", "num.ct.gt.1k", "perc.ct.gt.1k", "num.ct.gt.10k", "perc.ct.gt.10k", "num.ct.gt.100k", "perc.ct.gt.100k", "num.ct.gt.1m", "perc.ct.gt.1m", "num.ct.gt.10m", "perc.ct.gt.10m", "ct.mean.size", "ct.median.size", "ct.n50", "ct.l50", "ct.percA", "ct.percC", "ct.percG", "ct.percT", "ct.percN", "ct.perc.nonATGCN", "X")
cont<-cont[,-c(58,59)]

#get rid of non-reptiles
cont<-subset(cont, assembly != "homSap" & assembly != "danRer" & assembly != "musMus" & assembly != "xenTro")

#fix two errors
cont$assembly[cont$assembly=="ansCtg"]="ansCyg"
cont$assembly[cont$assembly=="chlUnd"]="chlMac"


#read in metadata
meta1<-read.table("../metadata//summary.out", header=T, sep="\t", stringsAsFactors=F)
meta2<-read.table("../metadata/assembly_result_04232015.txt", header=T, sep="\t", stringsAsFactors=F)
names(meta2)=c("genbank.acc", "genbank.id", "refseq.acc", "refseq.id", "drop")
meta2<-meta2[,-5]
meta3<-read.table("../metadata/genomes_euks_cleaned.txt", header=T, sep="\t", stringsAsFactors=F)
names(meta3)<-c("organism.name", "bioproject", "group", "subgroup", "size.mb", "gc", "assembly", "scaffold", "genes", "proteins", "init.date", "mod.date", "level")

library(XML)
library (plyr)
meta4<-xmlToList("../metadata/assembly_result_04232015.xml", simplify=T)
meta4.cat<-t(meta4)
meta4.cat<-meta4.cat[,c(12,6,3,35,29,31,15)]
meta4.df<-as.data.frame(matrix(unlist(meta4.cat), nrow=77, byrow=F), stringsAsFactors=F)
names(meta4.df)<-c("species", "asm.name", "acc", "submitter", "orig.date", "mod.date", "status")

#merge
#first, make a acc "database"
acc.key<-meta2
acc.key$acc<-acc.key$genbank.acc
acc.key=rbind(acc.key, acc.key)
acc.key$acc[78:154]=acc.key$refseq.acc[78:154]
acc.key=subset(acc.key, acc!="N/A")
genome.list<-read.table("../metadata//genome_versions_used.txt", sep=" ", header=F, stringsAsFactors=F)
genome.list=genome.list[,c(2,3)]
names(genome.list)=c("shortname", "acc")
acc.key<-merge(acc.key, genome.list, all=T)
acc.key=subset(acc.key, !is.na(shortname))
acc.key=rbind(acc.key, acc.key)
acc.key$acc[1:72]=acc.key$genbank.acc[1:72]
acc.key$acc[73:144]=acc.key$refseq.acc[73:144]
acc.key=subset(acc.key, !is.na(acc))

meta.merge<-merge(meta4.df, acc.key, by.x="acc", by.y="acc", all=T)
meta.merge<-subset(meta.merge, !is.na(shortname))
meta.merge<-subset(meta.merge, !is.na(species))

meta3$assembly=sub("\\s+$", "", meta3$assembly, perl=T) #clean up whitespace
meta.merge<-merge(meta.merge, meta3, by.x="asm.name", by.y="assembly", all.x=T, all.y=F)

makeShort<-function(x) {
  n<-strsplit(x, split="\\s+")
  name<-paste0(substr(tolower(n[[1]][1]), 1,3), substr(toupper(n[[1]][2]),1,1), substr(tolower(n[[1]][2]),2,3))  
  return(name)
}
meta1<-subset(meta1, !is.na(species))
meta1$shortname=apply(meta1[,c("species"),drop=F], 1, makeShort)

meta.merge<-merge(meta.merge, meta1, by="shortname", all=T)

#final merge
cont.all<-merge(cont, meta.merge, all=T, by.x="assembly", by.y="shortname")
cont.all$asm.name.y[is.na(cont.all$asm.name.y)]="N/A"

#clean up duplicates, remove hooded crow which I didn't include
cont.all<-cont.all[cont.all$asm.name.y != "AquilaChrysaetos1",]
cont.all<-subset(cont.all, asm.name.y != "SMACv1.0")
cont.all<-subset(cont.all, asm.name.y != "NB1.0")
cont.all<-subset(cont.all, assembly!="corCor")
cont.all<-subset(cont.all, select=-species.y)
cont.all<-subset(cont.all, select=-asm.name.y)
cont.all<-subset(cont.all, select=-c(mod.date.y, init.date))
names(cont.all)[c(58,60,63)]=c("asm.name", "species", "mod.date")

#add info on new genomes
our.genomes<-data.frame(assembly=cont.all$assembly[is.na(cont.all$asm.name)],submitter="Harvard",orig.date=c("2015/01/04"),mod.date=c("2015/01/04"),status="Scaffold",subgroup="Birds",method="ALLPATHS-LG",coverage="70x",seq.tech="Illumina HiSeq", stringsAsFactors=F)
our.genomes[13,]=c("strCam", "BGI", "2014/06/06 00:00", "2015/12/07", "Chromosome", "Birds", "SOAPdenovo v. 1.6","85x","Illumina HiSeq")
our.genomes[14,]=c("uriLom", "Harvard", "2015/04/20", "2015/04/20", "Scaffold", "Birds", "Platanus","80x","Illumina HiSeq")
our.genomes[8,]=c("haeMex", "Harvard", "2014/04/25", "2014/04/25", "Scaffold", "Birds", "ALLPATHS-LG","80x","Illumina HiSeq")

cont.all[is.na(cont.all$asm.name),colnames(our.genomes)]=our.genomes

#clean up submitter
cont.all$submitter[cont.all$submitter=="WUGSC" | cont.all$submitter=="The Genome Institute - Washington University School of Medicine" | cont.all$submitter=="The Genome Institute, Washington University at St. Louis" | cont.all$submitter=="Broad Institute" | cont.all$submitter == "Washington University Genome Sequencing Center" | cont.all$submitter=="BCM-HGSC"] = "Genome Center"
cont.all$submitter[cont.all$submitter=="The Consortium for Comparative Genomics, UC Denver"] = "Other Independent"
cont.all$submitter[grepl("consortium", cont.all$submitter, ignore.case=T)] = "Consortium"
cont.all$submitter[grepl("BGI", cont.all$submitter, ignore.case=T)] = "BGI"
cont.all$submitter[cont.all$submitter=="Beijing Genomics Institute"] = "BGI"
cont.all$submitter[cont.all$submitter=="The International Crocodilian Genomes Working Group" | cont.all$submitter=="International Crocodilian Genomes Working Group"] = "Consortium"
cont.all$submitter[cont.all$submitter != "Harvard" & cont.all$submitter != "BGI" & cont.all$submitter != "Genome Center" & cont.all$submitter != "Consortium"] = "Other Independent"

#clean up method
cont.all$method[grepl("SOAPdenovo", cont.all$method, ignore.case=T)]="SOAPdenovo"
cont.all$method[grepl("ALLPATHS", cont.all$method, ignore.case=T)]="ALLPATHS-LG"
cont.all$method[cont.all$method=="Soap deNovo v. March 2012"]="SOAPdenovo"
cont.all$method[grepl("Celera", cont.all$method, ignore.case=T)]="Celera"
cont.all$method[grepl("CLC", cont.all$method, ignore.case=T)]="CLC"

#meta data and contiguity metrics are now in cont.all

#make plot for lab meeting
cont.all$plotcolor="forestgreen"
cont.all$plotcolor[cont.all$method=="SOAPdenovo"] = "blue"
cont.all$plotcolor[cont.all$method=="ALLPATHS-LG"]="red"
cont.all$plotcolor[cont.all$method=="Platanus"] = "orange"
cont.all$plotcolor[cont.all$method=="CLC"] = "orange"
cont.all$plotcolor[cont.all$method=="Ray software v. 3"] = "orange"
cont.all$plotsym = 15
cont.all$plotsym[cont.all$subgroup != "Birds"] = 17

par(mar=c(5,8,2,2))
with(cont.all[cont.all$assembly != "haeMex" & cont.all$assembly != "uriLom" & cont.all$assembly != "oceLeu" & cont.all$subgroup == "Birds",], plot(scaf.n50, ct.n50, xlab="Scaffold N50", ylab="", log="xy", las=1, col=plotcolor, pch=16, cex=2, cex.axis=2, cex.lab=2))
mtext(2, text="Contig N50", 6.5, cex=2)
with(cont.all[cont.all$assembly != "haeMex" & cont.all$assembly != "uriLom" & cont.all$assembly != "oceLeu" & cont.all$subgroup == "Birds" & cont.all$submitter=="Harvard" & cont.all$assembly!="aptHaa" & cont.all$assembly!="aptRow" & cont.all$assembly!="aptOwe" & cont.all$assembly != "notPer" & cont.all$assembly != "droNov" & cont.all$assembly != "casCas",], text(scaf.n50, ct.n50, labels=assembly, pos=4, offset=0.4, cex=1.2))
with(cont.all[cont.all$assembly=="aptHaa" | cont.all$assembly=="aptOwe" | cont.all$assembly== "notPer" | cont.all$assembly == "droNov",], text(scaf.n50, ct.n50, labels=assembly, pos=2, offset=0.4, cex=1.2))
with(cont.all[cont.all$assembly=="aptRow",], text(scaf.n50, ct.n50, labels=assembly, adj=c(0,1.25), cex=1.2))
with(cont.all[cont.all$assembly=="casCas",], text(scaf.n50, ct.n50, labels=assembly, adj=c(.3,-.6), cex=1.2))
with(cont.all[cont.all$assembly=="strCam",], text(scaf.n50, ct.n50, labels=assembly, adj=c(.3,-.6), cex=1.2))
with(cont.all[cont.all$assembly=="tinGut",], text(scaf.n50, ct.n50, labels=assembly, adj=c(.3,-.6), cex=1.2))
with(cont.all[cont.all$assembly=="galGal",], text(scaf.n50, ct.n50, labels=assembly, adj=c(.8,-.6), cex=1.2))
with(cont.all[cont.all$assembly=="melGal",], text(scaf.n50, ct.n50, labels=assembly, adj=c(.8,-.6), cex=1.2))
with(cont.all[cont.all$assembly=="ficAlb",], text(scaf.n50, ct.n50, labels=assembly, adj=c(.8,-.6), cex=1.2))
with(cont.all[cont.all$assembly=="taeGut",], text(scaf.n50, ct.n50, labels=assembly, adj=c(.8,-.6), cex=1.2))


with(cont.all[cont.all$assembly=="pseHum",], text(scaf.n50, ct.n50, labels=assembly, adj=c(.8,-.6), cex=1.2))
legend("topleft", legend=c("ALLPATHS-LG", "SOAPdenovo", "other short-read", "Sanger"), col=c("red", "blue", "orange", "forestgreen"), pch=16)

##BUSCO
busco<-read.table("../completeness/busco_v1/all_busco_sum.txt", header=T)
all.genomes<-read.table("assembly_summary_genbank.txt", header=T, sep="\t", comment.char="", quote="", colClasses="character")
genomes<-read.table("bird_assemblies_at_genbank.txt", header=T, sep="\t", colClasses="character")
names(genomes)[1]="accession"
names(all.genomes)[1]="accession"
submitter=all.genomes[,c("accession", "submitter")]
genome.key=genomes[,c("accession", "shortname", "organism_name", "asm_name")]
genome.key=merge(genome.key, submitter, all.x=T)
genome.key=genome.key[genome.key$shortname != "strCam",]
genome.key=rbind(genome.key, data.frame(accession=c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"), 
                                        shortname=c("rheAme", "strCam", "casCas", "haeMex", "droNov", "aptOwe", "aptHaa", "cryCin", "oceLeu"), 
                                        organism_name=c("Rhea americana", "Struthio camelus", "Casuarius casuarius", "Haemorhous mexicanus", "Dromaius novaehollandiae", "Apteryx owenii", "Apteryx haastii", "Crypturellus cinnamomeus", "Oceanodroma leucorhoa"),
                                        asm_name=c("rheAme1", "strCamOM", "casCas1", "haeMex1", "droNov1", "aptOwe1", "aptHaa1", "cryCin1", "oceLeu1"),
                                        submitter=c("Team Ratite", "Bachtrog", "Team Ratite", "Team HouseFinch", "Team Ratite", "Team Ratite", "Team Ratite", "Team Ratite", "Team Petrel")))
#clean up submitters
genome.key$submitter[genome.key$submitter=="Beijing Genomics Institute"]="BGI"
genome.key$submitter[genome.key$submitter=="The Genome Institute - Washington University School of Medicine"]="WashU"
genome.key$submitter[genome.key$submitter=="The Genome Institute,  Washington University at St. Louis"]="WashU"
genome.key$submitter[genome.key$submitter=="Beijing Genomics Institute (BGI)-Shenzhen"]="BGI"
genome.key$submitter[genome.key$submitter=="Washington University Genome Sequencing Center"]="WashU"

busco=merge(busco, genome.key, by.x="species", by.y="shortname", all.x=T, all.y=F)

busco.sum<-as.data.frame(table(busco$sp, busco$len < 50))
busco.sum=busco.sum[busco.sum$Var2==TRUE,]
busco.sum$missing=busco.sum$Freq/2882
busco.sum$missing[busco.sum$Var1 == "ficAlb"]=busco.sum$Freq[busco.sum$Var1 == "ficAlb"]/332
busco.sum$missing.plot = busco.sum$missing
busco.sum$missing.plot[busco.sum$missing>0.35] = 0.35

busco$hit = cut(busco$len, breaks=c(-1, 50, 70, 101), labels=c("missing", "patial", "present"))
busco.sum2<-as.data.frame(table(busco$sp, busco$hit)/2882)
busco.sum2$Freq[busco.sum2$Var1=="ficAlb"]=(busco.sum2$Freq[busco.sum2$Var1=="ficAlb"]*2882)/332
busco.barplot<-unstack(busco.sum2, Freq~Var2)
busco.barplot$species=busco.sum2$Var1[busco.sum2$Var2=="missing"]


plot(busco.sum$missing.plot[order(busco.sum$missing)], xaxt="n", yaxt="n", xlab="", ylab="Fraction Missing", ylim=c(0,0.35))
axis(1, labels=c(as.character(busco.sum$Var1[order(busco.sum$missing)])), at=c(1:length(busco.sum$Var1)), las=2)
axis(2, at=c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35), las=1)

barplot(t(busco.barplot[order(busco.barplot$missing),1:3]), names=busco.barplot$species[order(busco.barplot$missing)], las=2, legend.text=c("missing", "partial", "present"), args.legend=c(x=8.4,y=1), ylab="Fraction")

#plot
new.ratites<-data.frame(Var1 = c("aptHaa", "aptOwe", "aptRow", "casCas", "cryCin", "droNov", "eudEle", "notPer", "rheAme", "rhePen", "strCam"), organism.name = c("Apteryx haastii","Apteryx owenii","Apteryx rowi","Casuarius casuarius","Crypturellus cinnamomeus", "Dromaius novaehollandiae", "Eudromia elegans", "Nothoprocta perdicaria", "Rhea americana", "Rhea pennata", "Struthio camelus"))

name.key = subset(cont.all, !is.na("organism.name"), select=c("assembly", "organism.name"))
names(name.key)[1]="Var1"
name.key = rbind(name.key, new.ratites)

par(mar=c(12,6,0.5,0.5))
busco.plot=merge(busco.sum,name.key) %>% filter(!is.na(organism.name))
busco.plot$missing[busco.plot$missing > 0.30] = 0.30
busco.plot$plotcol="skyblue"
busco.plot$plotcol[busco.plot$Var1 %in% c("droNov", "aptHaa", "aptRow", "aptOwe", "casCas", "rheAme", "rhePen", "eudEle", "notPer", "cryCin")]="red"
busco.plot$plotcol[busco.plot$Var1 %in% c("melGal", "galGal", "ficAlb", "taeGut")]="purple"
plot(busco.plot$missing[order(busco.plot$missing, decreasing=F)], xaxt="n", xlab="", ylab="", col="gray20", type="h", las=2, ylim=c(0,0.3), cex.axis=2, lwd=1)
points(busco.plot$missing[order(busco.plot$missing, decreasing=F)], pch=ifelse(busco.plot$missing[order(busco.plot$missing, decreasing=F)]>=0.3, 17, 16), type="p", col=busco.plot$plotcol[order(busco.plot$missing, decreasing=F)], cex=1.2)
axis(1, labels=c(as.character(busco.plot$organism.name[order(busco.plot$missing, decreasing=F)])), at=c(1:length(busco.plot$Var1)), las=2, cex.axis=1)

legend(x="topleft", legend=c("New Palaeognaths", "Short Read", "\"Finished\" Genomes"), col=c("red", "skyblue", "purple", "black"), pch=16, cex=2, bty="n")
mtext(2, text="Fraction Highly Conserved Genes Absent", line=4.6, cex=2)


