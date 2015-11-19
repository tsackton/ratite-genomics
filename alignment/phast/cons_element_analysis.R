#Analyze conserved elements

ce1<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/most_conserved_final.tree1.bed", header=F)
ce2<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/most_conserved_final.tree2.bed", header=F)
lowe<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/final_beds/lowe_cnees.bed", header=F)

#add lengths to each
ce1$length<-ce1$V3-ce1$V2
ce2$length<-ce2$V3-ce2$V2
lowe$length<-lowe$V3-lowe$V2

#plot length distributions
plot(density(ce2$length[ce2$length<1000]), xlim=c(0,300), col="red", lwd=2, main="Conserved element length distribution", xlab="Length", bty="n")
lines(density(lowe$length[lowe$length<1000]), lwd=2, lty="dashed")
legend("topright", legend=c("Sackton et al", "Lowe et al"), lty=c("solid", "dashed"), col=c("red", "black"), lwd=2, bty="n")
abline(v=median(ce2$length), col="red", lty="dotted", lwd=2)
abline(v=median(lowe$length), col="black", lty="dotted", lwd=2)

#length distibutions of elements > 50 bp
plot(density(ce2$length[ce2$length<1000 & ce2$length > 75]), xlim=c(0,300), col="red", lwd=2, main="Conserved element length distribution (75bp cutoff)", xlab="Length", bty="n")
lines(density(lowe$length[lowe$length<1000 & lowe$length > 75]), lwd=2, lty="dashed")
legend("topright", legend=c("Sackton et al", "Lowe et al"), lty=c("solid", "dashed"), col=c("red", "black"), lwd=2, bty="n")
abline(v=median(ce2$length[ce2$length<1000 & ce2$length > 75]), col="red", lty="dotted", lwd=2)
abline(v=median(lowe$length[lowe$length<1000 & lowe$length > 75]), col="black", lty="dotted", lwd=2)

plot(hexbin(ce2$length, ce2$V5))
