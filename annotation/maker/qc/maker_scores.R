species=c("aptHaa", "aptRow", "aptOwe", "droNov", "casCas", "rheAme", "rhePen", "cryCin", "notPer", "eudEle")

maker.scores<-list()
maker.qc<-list()

for (sp in species) {
  aed.in<-read.table(paste0("~/Projects/birds/ratite_compgen/ratite-genomics/annotation/maker/qc/", sp, ".transcripts.info"), sep=" ");
  #clean up
  aed.in$gene = sub(">(\\w+)-R\\w", "\\1", aed.in$V1, perl=T)
  aed.in$isoform = sub(">\\w+-(\\w+)", "\\1", aed.in$V1, perl=T)
  aed.in$aed = sub("AED:", "", aed.in$V4)
  aed.in$eaed = sub("eAED:", "", aed.in$V5)
  aed.in$qi=sub("QI:", "", aed.in$V6)
  qistring = strsplit(aed.in$qi, "|", fixed=T)
  qistring2 = do.call("rbind", qistring)
  aed.in[,12:20]<-qistring2
  names(aed.in)[12:20]<-c("len.5utr", "ss.est", "exon.est", "exon.prot", "ss.snap", "exon.snap", "exon.num", "len.3utr", "len.prot")
  aed.fix<-aed.in$eaed
  maker.qc[[sp]]=aed.in[,7:20]
  maker.scores[[sp]]=data.frame(AED=sort(as.numeric(as.character((aed.fix)))), cum.frac=seq(1,length(aed.fix))/length(aed.fix), stringsAsFactors=F)
  
}

#plot
plot(x=1, y=1, xlim=c(0,1), ylim=c(0,1), xlab="Annotation Edit Distance", ylab="Cumulative Fraction", type="n", bty="l", cex.axis=1, las=1)
ratite.colors=data.frame(species=species, col=c("saddlebrown", "saddlebrown", "saddlebrown", "royalblue2", "red3", "seagreen", "seagreen", "plum2", "plum2", "plum2"), stringsAsFactors=F)
for (sp in species) {
  lines(maker.scores[[sp]]$AED, maker.scores[[sp]]$cum.frac, col=ratite.colors$col[ratite.colors$sp==sp], lwd=2)
}
legend("bottomright", cex=2, legend=c("Kiwi", "Emu", "Cassowary", "Rhea", "Tinamou"), col=c("saddlebrown", "royalblue2", "red3", "seagreen", "plum2"), lwd=2)
abline(v=0.5, lty="dashed", lwd=2)

#number of genes
lapply(maker.qc, function(x) length(unique(x$gene)))

maker_clean <- do.call("rbind", maker.qc) %>% add_rownames %>% tbl_df %>%
  mutate(species = sub("\\.\\d+", "", rowname)) %>% select(species, gene:len.prot)

gene_sum <- maker_clean %>% group_by(species, gene) %>% summarize(isoform_num = n(), aed = max(as.numeric(eaed)))
gene_sum <- gene_sum %>% mutate(quality_bin = cut(aed, breaks=c(0,0.25,0.5,0.75,1), labels = F))
par(mar=c(6,4,4,1))
barplot(table(gene_sum$quality_bin, gene_sum$species), ylab="Number of genes")
legend(3,-2000, legend=c("AED 0-0.25", "AED 0.26-0.50", "AED 0.51-0.75", "AED 0.76-1"), col=c("gray20", "gray40", "gray60", "gray80"), pch=15, ncol=4, xpd=TRUE)
