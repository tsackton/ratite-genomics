species=c("aptHaa", "aptRow", "aptOwe", "droNov", "casCas", "rheAme", "rhePen", "cryCin", "notPer", "eudEle")

maker.scores<-list()

for (sp in species) {
  aed.in<-read.table(paste0("~/Projects/birds/ratite_compgen/ratite-genomics/annotation/maker/qc/", sp, ".transcripts"), sep=" ");
  aed.fix<-sub("AED:", "", aed.in$V4)
  maker.scores[[sp]]=data.frame(AED=sort(as.numeric(as.character((aed.fix)))), cum.frac=seq(1,length(aed.fix))/length(aed.fix), stringsAsFactors=F)
}

#plot
plot(x=1, y=1, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", type="n", bty="l", cex.axis=2, las=1)
ratite.colors=data.frame(species=species, col=c("saddlebrown", "saddlebrown", "saddlebrown", "royalblue2", "red3", "seagreen", "seagreen", "plum2", "plum2", "plum2"), stringsAsFactors=F)
for (sp in species) {
  lines(maker.scores[[sp]]$AED, maker.scores[[sp]]$cum.frac, col=ratite.colors$col[ratite.colors$sp==sp], lwd=4)
}
legend("bottomright", cex=2, legend=c("Kiwi", "Emu", "Cassowary", "Rhea", "Tinamou"), col=c("saddlebrown", "royalblue2", "red3", "seagreen", "plum2"), lwd=4)
