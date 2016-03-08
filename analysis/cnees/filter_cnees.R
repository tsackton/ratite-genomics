#filtering for CNEE candidates

conv<-read.table("~/Projects/birds/ratite_compgen/ratite-genomics/alignment/moa/raw_convergence_results.tsv.gz", header=T, stringsAsFactors=F)

#subset CNEE list to exclude low coverage or duplicated cases
#these filters are kind of arbitrary... (75% coverage of both birds and tinamous, <= 3 duplicated)
conv.good<-subset(conv, AvesPres >= 27 & tinPres >= 3 & Freq <= 3)

#now make a separate data frame for each of kiwi, moa, casaur, strcam, rhea, also tinamou accel and basalPaleo accel

kiwi<-spread(subset(conv, select=c("name", "model", "kiwi.br")), model, kiwi.br)
ostrich<-spread(subset(conv, select=c("name", "model", "strcam.br")), model, strcam.br)
rhea<-spread(subset(conv, select=c("name", "model", "rhea.br")), model, rhea.br)
casuar<-spread(subset(conv, select=c("name", "model", "casuar.br")), model, casuar.br)
moa<-spread(subset(conv, select=c("name", "model", "moa.br")), model, moa.br)

kiwi$new = apply(kiwi[,c(2,3,4)], 1, max)
ostrich$new = apply(ostrich[,c(2,3,4)], 1, max)
rhea$new = apply(rhea[,c(2,3,4)], 1, max)
casuar$new = apply(casuar[,c(2,3,4)], 1, max)
moa$moa.c = apply(moa[,c(2,3,4)], 1, max)
moa=subset(moa, select=c("name", "moa.c"))

kiwi$kiwi.c<-kiwi$original*kiwi$new
kiwi$kiwi.df<-abs(kiwi$original - kiwi$new)

ostrich$ostrich.c<-ostrich$original*ostrich$new
ostrich$ostrich.df<-abs(ostrich$original - ostrich$new)

rhea$rhea.c<-rhea$original*rhea$new
rhea$rhea.df<-abs(rhea$original - rhea$new)

casuar$casuar.c<-casuar$original*casuar$new
casuar$casuar.df<-abs(casuar$original - casuar$new)

filter<-cbind(kiwi[,c("name", "kiwi.c", "kiwi.df")], ostrich[,c("ostrich.c", "ostrich.df")], rhea[,c("rhea.c", "rhea.df")], casuar[,c("casuar.c", "casuar.df")])

filter$moa.c=moa$moa.c

filter$max.df = apply(filter[,c(3,5,7,9)], 1, max)
filter$moa.c[filter$max.df >0] = 0

filter$clade.ct = filter$kiwi.c + filter$ostrich.c + filter$moa.c + filter$rhea.c + filter$casuar.c

basal<-spread(subset(conv, select=c("name", "model", "basalPaleo")), model, basalPaleo)
basal$b.min = apply(basal[,c(2,3,4)], 1, min)

tinamou<-spread(subset(conv, select=c("name", "model", "tinamou")), model, tinamou)
tinamou$t.min = apply(tinamou[,c(2,3,4)], 1, min)

accel<-cbind(tinamou[,c("name", "t.min")], basal[,c("b.min")])
accel$palaeo.test = apply(accel[,c(2,3)], 1, min)

filter$accel.palaeo = accel$palaeo.test

filter.best = subset(filter, accel.palaeo > 0.25)

cnee.info = unique(conv[,c("name", "best_ens", "best_ncbi")])
filter.final=merge(filter.best, cnee.info, all=F)
