#code for plotting CNEEs

#first read in file
#change the path to wherever you have the file stored
cnee<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/phast/ratite_data_ver1.tsv", header=T, sep="\t")

#add some columns
cnee$ratite_accel = cnee$e_and_r + cnee$r_only
cnee$emu_accel = cnee$e_and_r + cnee$e_only

#get fractions
cnee$total_frac = cnee$total.accel / cnee$total
cnee$ratite_frac = cnee$ratite_accel / cnee$total
cnee$emu_frac = cnee$emu_accel / cnee$total

#subset to only keep genes with cnees
#there are some genes with total = 0 meaning we didn't find 
#any cnees > 50 bp near that gene
cnee.2 = subset(cnee, total > 0)

#define a cutoff for the plot; everything with a frac > cutoff will be a different color
cutoff = 0.2

plot(total ~ ratite_accel, log="y", data=cnee.2, las=1, pch=16, cex=0.5, col=ifelse(ratite_frac > cutoff, "red", "black"), xlab="# CNEEs accelerated in ratites", ylab="Total # of CNEEs")

#add lines to show the % cutoffs
abline(a=0, b=3, untf=T, lty="dashed", lwd=2, col="blue") # 33% of CNEEs accelerated
abline(a=0, b=5, untf=T, lty="dashed", lwd=2, col="blue") # 20%
abline(a=0, b=10, untf=T, lty="dashed", lwd=2, col="blue") # 10%

#add text
text(x=55, y=210, labels=("20% accel"))
text(x=60, y=450, labels=c("10% accel"))
text(x=50, y=80, labels=c("33% accel"))

