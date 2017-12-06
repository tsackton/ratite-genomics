setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/protein_coding/paml_branch/")
library(data.table)
library(qvalue)
library(tidyverse)
library(clusterProfiler)
library(org.Gg.eg.db)

##SETUP##

#read in hog_info key
#read in ancrec parsing key for species tree info

ancrec.parsed<-fread("gunzip -c ../paml_ancrec/paml_M0_parsed.txt.gz")
ancrec.treekey<-ancrec.parsed[,c("hog", "treenum", "species_tree"), with=FALSE]
hog_info<-read.table("../all_hog_info.tsv", sep="\t", header=T)

hog_info$has_species_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree]
hog_info$has_gene_tree = hog_info$hog %in% ancrec.treekey$hog[ancrec.treekey$species_tree == F]

#load hog <-> chicken info
hog_to_gene <- read.table("../HOG_final_alignment_seqids", header=T, stringsAsFactors =F)
hog_to_gene <- hog_to_gene %>% tbl_df %>%
  mutate(hog = as.integer(sub("HOG2_", "", HOG, fixed=T))) %>% 
  filter(Taxon == "galGal") %>%
  select(-HOG, -Fasta_seqid, -Taxon, -Transcript) %>%
  group_by(hog) %>%
  mutate(gene = paste0(Gene, sep="", collapse=";")) %>%
  select(-Gene)

#define clades -- fixed (vocal learners, palaeognaths, tinamous, ratites)
all_clades<-names(hog_info)[2:40]
ratite_clades<-c("aptHaa", "aptOwe", "aptRow", "rheAme", "rhePen", "strCam", "casCas", "droNov")
vl_clades<-c('calAnn', 'corBra', 'serCan', 'geoFor', 'melUnd', 'pseHum', 'taeGut', 'ficAlb')
tin_clades<-c("tinGut", "cryCin", "notPer", "eudEle")
palaeo_clades<-c(ratite_clades, tin_clades)
non_ratite_clades<-all_clades[!(all_clades %in% ratite_clades)]
non_ratite_palaeo<-palaeo_clades[!(palaeo_clades) %in% ratite_clades]

##ANALYSIS - FUNCTIONS##

#setup - load and clean data
prep_data <- function(file, ancrec.treekey, hog_info, check_missing = F) {
  #load data -- dn
  read_line = paste0("gunzip -c ", file)
  dn<-fread(read_line, header=F, sep=",")
  names(dn)<-c("hog", "tree", "parent.node", "desc.node", "branch.id", "dn", "ratite", "vl")
  
  #add hog_info
  dn$tree = sub("tree", "", dn$tree, fixed=T)
  dn$tree = as.integer(dn$tree)
  dn<-merge(dn, ancrec.treekey, by.x=c("hog", "tree"), by.y=c("hog", "treenum"), all.x=T, all.y=F)
  dn<-merge(dn, hog_info, by.x="hog", by.y="hog", all=T)
  #get missing runs
  if (check_missing) {
    check_for_missing<-unique(subset(dn, is.na(dn), select=c("hog", "has_species_tree", "has_gene_tree")))
    write.table(check_for_missing$hog, file="branch_hogs_torerun.txt", quote=F, sep="", row.names = F, col.names =  T)
    print(nrow(check_for_missing))
   }
  return(dn)
}

subset_clean_data <- function(DF, missing_cutoff = 2, dup_cutoff = 0, use_sptree = TRUE) {
  dn.clean = subset(DF, dup_ct <= dup_cutoff & missing_ct <= missing_cutoff  & species_tree == use_sptree, select=c("hog", "parent.node", "desc.node", "branch.id", "dn"))
  dn.clean$ratite=sapply(strsplit(dn.clean$desc.node, "-"), function(x) sum(x %in% ratite_clades)/length(x) == 1)
  dn.clean$vl=sapply(strsplit(dn.clean$desc.node, "-"), function(x) sum(x %in% vl_clades)/length(x) == 1)
  dn.clean$nrpalaeo=sapply(strsplit(dn.clean$desc.node, "-"), function(x) sum(x %in% non_ratite_palaeo)/length(x) == 1)
  return(dn.clean)
}

#function to do vector projection
proj_vect <- function(genevec, sptree) {
  as.matrix(genevec) - sptree %*% t(sptree) %*% as.matrix(genevec)
}
#actually compute normalized stat
normalize_branch_stat <- function(DF, filter=TRUE) {
  dn.clean<-as.data.table(DF)
  #make unit vector
  dn.clean[,dn.length.bygene:=sqrt(sum(dn^2)), by=list(hog)]
  dn.clean$dn.unit.bygene=dn.clean$dn/dn.clean$dn.length.bygene
  #make average tree (average of all branches a tree appears on)
  dn.clean[,dn.average.tree:=mean(dn.unit.bygene, na.rm=T), by=list(branch.id)]
  #convert to unit vector (this will be different for each species tree configuration)
  dn.clean[,dn.unit.sptree:=dn.average.tree/sqrt(sum(dn.average.tree^2)), by=.(hog)]
  dn.clean[,dn.norm := proj_vect(dn.unit.bygene, dn.unit.sptree), by=.(hog)]
  if (filter) {
    branch_freqs<-as.data.frame(table(dn.clean$branch.id))
    dn.clean<-dn.clean[dn.clean$branch.id %in% branch_freqs$Var1[branch_freqs$Freq >= 500],]
  }
  return(dn.clean)
}

compute_pval <- function (x, groupvar) {
  if (inherits(try(ans<-wilcox.test(x ~ groupvar, conf.int=TRUE),silent=TRUE),"try-error"))
    return(NA_real_)
  else
    return(list(pval=ans$p.value, est=ans$estimate))
}

#qc checks of distributions
make_qc_plots <- function(DF, file) {
  pdf(file=file)
  hist(DF$dn.norm, breaks=100)
  branch_freqs<-as.data.frame(table(dn.clean$branch.id))
  #normalization checks
  high_freq_branches = branch_freqs$Var1[branch_freqs$Freq > 5000]
  branch_key = unique(dn.clean[,c("branch.id","ratite","vl"), with=F])
  plot_branch_dists = data.frame(dn.norm=dn.clean$dn.norm[dn.clean$branch.id %in% high_freq_branches], branch.id=dn.clean$branch.id[dn.clean$branch.id %in% high_freq_branches])
  plot_branch_dists = merge(plot_branch_dists, branch_key, by.x="branch.id", by.y="branch.id")
  plot_branch_dists$color = "black"
  plot_branch_dists$color[plot_branch_dists$ratite]="red"
  plot_branch_dists$color[plot_branch_dists$vl]="blue"
  plot_branch_dists_colors = unique(plot_branch_dists[,c("branch.id", "color")])
  boxplot(plot_branch_dists$dn.norm2 ~ plot_branch_dists$branch.id, outline=F, col=plot_branch_dists_colors$color)
  dev.off()
}

do_perms<-function(DF, nreps=1000, load="", write="") {
  #permutations
  if (load != "") {
    dn.perm<-fread(load)
    return(dn.perm)
  }
  
  dn.perm=DF[,.(hog,desc.node,dn.norm)]
  dn.perm.desc_nodes=data.frame(desc.node=unique(dn.perm$desc.node),test_clade=F, stringsAsFactors = F)
  nreps=nreps
  for (i in 1:nreps) {
    clades_okay = F
    while (!clades_okay) {
      test_clades = sample(all_clades, 8)
      clades_okay = sum(test_clades %in% ratite_clades) < 2 & sum(test_clades %in% vl_clades) < 2 
    }
    dn.perm.desc_nodes$test_clade=sapply(strsplit(dn.perm.desc_nodes$desc.node, "-"), function(x) sum(x %in%test_clades)/length(x) == 1)
    dn.perm$test_clade = dn.perm$desc.node %in% dn.perm.desc_nodes$desc.node[dn.perm.desc_nodes$test_clade]
    newcol=paste0("rep",i)
    dn.perm[,(newcol):= { if (inherits(try(ans<-wilcox.test(dn.norm ~ test_clade)$p.value,silent=TRUE),"try-error"))
      NA_real_
      else
        ans
    }, by=list(hog)]
  }
  
  dn.perm.pval=unique(dn.perm[,c(1,5:length(dn.perm)), with=F])
  if (write != "") {
    fwrite(dn.perm.pval, file=write)
  }
  return(dn.perm.pval)
}

get_dir <- function(x, groupby) {
  cor(x=x, y=as.numeric(groupby), use="na.or.complete")
}

get_coef <- function(x, groupby, subset) {
  group <- factor(groupby, levels=c("TRUE", "FALSE"))
  if (length(levels(droplevels(group))) == 2) {
    wilcox.test(x ~ group, conf.int=T)$estimate
  }
  else {
    return(NA_real_)
  }
}

## ANALYSIS STARTS HERE ###

dn<-prep_data(file="aa_parsed.csv.gz", ancrec.treekey = ancrec.treekey, hog_info = hog_info, check_missing = T)
dn.default<-subset_clean_data(dn)
#fix vl data, don't want to consider branches with descendant nodes that include the base of passerines/parrots as vocal learners
dn.default$vl[grepl("-melUnd", dn.default$desc.node, fixed=T)]=FALSE
dn.default$ratite[grepl("aptHaa-aptOwe-aptRow-casCas-droNov-rheAme-rhePen", dn.default$desc.node, fixed=T)]=FALSE
dn.default$ratite[grepl("aptHaa-aptOwe-aptRow-casCas-droNov", dn.default$desc.node, fixed=T)]=FALSE

#test out other convergences

dn.default$wb=FALSE
dn.default$wb[grepl("anaPla", dn.default$desc.node, fixed=T)]=TRUE
dn.default$wb[grepl("aptFor", dn.default$desc.node, fixed=T)]=TRUE
dn.default$wb[grepl("chaVoc", dn.default$desc.node, fixed=T)]=TRUE
dn.default$wb[grepl("egrGar", dn.default$desc.node, fixed=T)]=TRUE
dn.default$wb[grepl("nipNip", dn.default$desc.node, fixed=T)]=TRUE
dn.default$wb[grepl("pygAde", dn.default$desc.node, fixed=T)]=TRUE
dn.default$wb[grepl("balReg", dn.default$desc.node, fixed=T)]=TRUE
dn.default$wb[grepl("fulGla", dn.default$desc.node, fixed=T)]=TRUE
dn.default$wb[grepl("pseHum", dn.default$desc.node, fixed=T)]=FALSE
dn.default$wb[grepl("galGal-melGal-anaPla", dn.default$desc.node, fixed=T)]=FALSE
dn.default$wb[grepl("aptFor-pygAde-fulGla-egrGar-nipNip", dn.default$desc.node, fixed=T)]=FALSE

dn.default$tinamou=FALSE
dn.default$tinamou[grepl("notPer", dn.default$desc.node, fixed=T)]=TRUE
dn.default$tinamou[grepl("tinGut", dn.default$desc.node, fixed=T)]=TRUE
dn.default$tinamou[grepl("cryCin", dn.default$desc.node, fixed=T)]=TRUE
dn.default$tinamou[grepl("eudEle", dn.default$desc.node, fixed=T)]=TRUE
dn.default$tinamou[grepl("aptHaa", dn.default$desc.node, fixed=T)]=FALSE
dn.default$tinamou[grepl("cryCin-tinGut-eudEle-notPer", dn.default$desc.node, fixed=T)]=FALSE

dn.default<-normalize_branch_stat(dn.default)

dn.results<-dn.default[,.(ratite.coef = get_coef(dn.norm, ratite), 
                          nrpalaeo.coef = get_coef(dn.norm, nrpalaeo), 
                          ratite.pval = compute_pval(dn.norm, ratite), 
                          nrpalaeo.pval = compute_pval(dn.norm, nrpalaeo), 
                          wb.coef = get_coef(dn.norm, wb), 
                          vl.coef = get_coef(dn.norm, vl), 
                          wb.pval = compute_pval(dn.norm, wb), 
                          vl.pval = compute_pval(dn.norm, vl)), 
                       by=hog]
length(unique(dn.results$hog))

sum(dn.results$ratite.pval < 0.05 & dn.results$ratite.coef < 0) / sum(dn.results$ratite.coef < 0)
sum(dn.results$ratite.pval < 0.05 & dn.results$ratite.coef > 0) / sum(dn.results$ratite.coef > 0)

sum(dn.results$nrpalaeo.pval < 0.05 & dn.results$nrpalaeo.coef < 0, na.rm=T) / sum(dn.results$nrpalaeo.coef < 0, na.rm=T)
sum(dn.results$nrpalaeo.pval < 0.05 & dn.results$nrpalaeo.coef > 0, na.rm=T) / sum(dn.results$nrpalaeo.coef > 0, na.rm=T)

sum(dn.results$wb.pval < 0.05 & dn.results$wb.coef < 0) / sum(dn.results$wb.coef < 0)
sum(dn.results$wb.pval < 0.05 & dn.results$wb.coef > 0) / sum(dn.results$wb.coef > 0)

sum(dn.results$vl.pval < 0.05 & dn.results$vl.coef < 0, na.rm=T) / sum(dn.results$vl.coef < 0, na.rm=T)
sum(dn.results$vl.pval < 0.05 & dn.results$vl.coef > 0, na.rm=T) / sum(dn.results$vl.coef > 0, na.rm=T)

#simple perm
simple_perm <- function(perm, DF) {
  key<-DF %>% dplyr::select(desc.node, ratite) %>% distinct()
  key$perm = sample(key$ratite)
  dn.perm <- inner_join(key, dn.default) %>% as.data.table
  dn.perm<-dn.perm[,compute_pval(dn.norm, perm), by=hog] %>% dplyr::select(est, pval)
}

test_perm <- lapply(1:3, simple_perm, DF=dn.default)

#ratite, vocal learning, tinamou, waterbirds
#need to think of a better way of quantifying the signal in the data somehow

###TESTING

test_hog = 1293
test_df <- dn.default %>% filter(hog==test_hog)

test_df %>% arrange(dn.norm) %>% mutate(rank=rank(dn.norm), index=as.factor(ratite+nrpalaeo+tinamou)) %>% ggplot(aes(x=rank, y=dn.norm, col=index)) + geom_point() 

cor(test_df$dn.norm, test_df$ratite, method="sp")

mean(test_df$dn.norm[test_df$ratite]) - mean(test_df$dn.norm[!test_df$ratite])

library("boot")

boot_coef <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- rfit(formula, data=d)
  return(coef(fit)[2])
}

test_boot <- boot(data=test_df, statistic = boot_coef, R=1000, formula = dn.norm ~ ratite)

library("coin")

oneway_test(dn.norm ~ factor(, levels=c("TRUE", "FALSE")), data=test_df, distribution=approximate(B=9999))

library("Rfit")
fit<-rfit(dn.norm ~ ratite, data=test_df)

#below here needs editing


dn.pval<-dn.default[,.(ratite.p = compute_pval(dn.norm, ratite),ratite.dir = get_dir(dn.norm, ratite),wb.p=compute_pval(dn.norm, wb),wb.dir=get_dir(dn.norm, wb),vl.p=compute_pval(dn.norm, vl), vl.dir=get_dir(dn.norm,vl)), by=hog]
length(dn.pval$hog)
summary(qvalue(dn.pval$ratite.p))

dn.pval %>% filter(ratite.dir > 0) %>% ggplot(aes(x=ratite.p)) + geom_density(bins=50)
dn.pval %>% filter(ratite.dir > 0) %>% with(., summary(qvalue(ratite.p)))
dn.pval %>% filter(ratite.dir < 0) %>% ggplot(aes(x=ratite.p)) + geom_histogram(bins=50)
dn.pval %>% filter(ratite.dir < 0) %>% with(., summary(qvalue(ratite.p)))

dn.pval %>% filter(wb.dir > 0) %>% ggplot(aes(x=wb.p)) + geom_histogram(bins=50)
dn.pval %>% filter(wb.dir > 0) %>% with(., summary(qvalue(wb.p)))
dn.pval %>% filter(wb.dir < 0) %>% ggplot(aes(x=wb.p)) + geom_histogram(bins=50)
dn.pval %>% filter(wb.dir < 0) %>% with(., summary(qvalue(wb.p)))

dn.pval %>% filter(vl.dir > 0) %>% ggplot(aes(x=vl.p)) + geom_histogram(bins=50)
dn.pval %>% filter(vl.dir > 0) %>% with(., summary(qvalue(vl.p)))
dn.pval %>% filter(vl.dir < 0) %>% ggplot(aes(x=vl.p)) + geom_histogram(bins=50)
dn.pval %>% filter(vl.dir < 0) %>% with(., summary(qvalue(vl.p)))


#go enrich

accel_hogs <- dn.pval %>% filter(vl.dir > 0, vl.p < 0.01) %>% dplyr::select(hog) %>% left_join(hog_to_gene) %>% filter(!is.na(gene))
all_hogs <- dn.pval %>% dplyr::select(hog) %>% left_join(hog_to_gene) %>% filter(!is.na(gene))

enrichGO(accel_hogs$gene, 'org.Gg.eg.db', keytype = "SYMBOL", universe = all_hogs$gene, ont="MF")
enrichGO(accel_hogs$gene, 'org.Gg.eg.db', keytype = "SYMBOL", universe = all_hogs$gene, ont="BP")
enrichKEGG(bitr(accel_hogs$gene, fromType = "SYMBOL", toType="ENTREZID", OrgDb = "org.Gg.eg.db")$ENTREZID, organism="gga", keyType = "ncbi-geneid", universe = bitr(all_hogs$gene, fromType = "SYMBOL", toType="ENTREZID", OrgDb = "org.Gg.eg.db")$ENTREZID)

dn.perm.pval<-do_perms(dn.default, nreps=5)

#repeating subset and pval with different filtering
for (missing in c(0,2,4)) {
  for (usesp in c(TRUE, FALSE)) {
    print(missing)
    print(usesp)
    dn.default<-subset_clean_data(dn, missing_cutoff = missing, use_sptree = usesp)
    dn.default<-normalize_branch_stat(dn.default, filter=usesp)
    dn.pval<-dn.default[,.(ratite.p = compute_pval(dn.norm, ratite),vl.p=compute_pval(dn.norm, vl)), by=hog]
    print(length(dn.pval$hog))
    summary(qvalue(dn.pval$ratite.p))
    summary(qvalue(dn.pval$vl.p))
  }
}

#okay now let's look at dropping species and see how things change
#use only hogs with no missing data
dn.drop<-subset_clean_data(dn, missing_cutoff = 0)
#how many hogs?
dn.drop %>% distinct(hog) %>% summarize(count=n())
dn.drop <- normalize_branch_stat(dn.drop)

vl_tips_p = c("corBra", "ficAlb", "pseHum", "serCan", "taeGut", "geoFor")
vl_tips_np = c("melUnd", "calAnn")
vl_tips <- c("corBra", "ficAlb", "pseHum", "serCan", "taeGut", "geoFor", "melUnd", "calAnn")

get_new_target <- function(drop, nodekey, orig) {
  vl_tips <- c("corBra", "ficAlb", "pseHum", "serCan", "taeGut", "geoFor", "melUnd", "calAnn")
  tips_to_keep <- vl_tips[!(vl_tips %in% drop)]
  new_target<-logical(length(orig))
  for (selTip in tips_to_keep) {
    new_target[grepl(selTip, nodekey, fixed=T) & orig]=TRUE
  }
  return(new_target)
}

#make a bunch of new columns that each drop 1 vocal learner tips -- start there at least
dn.drop <- dn.drop %>% mutate(vl1p1 = get_new_target(c("corBra"), desc.node, vl)) %>%
  mutate(vl1p2 = get_new_target(c("ficAlb"), desc.node, vl)) %>%
  mutate(vl1p3 = get_new_target(c("pseHum"), desc.node, vl)) %>%
  mutate(vl1p4 = get_new_target(c("serCan"), desc.node, vl)) %>%
  mutate(vl1p5 = get_new_target(c("taeGut"), desc.node, vl)) %>%
  mutate(vl1p6 = get_new_target(c("geoFor"), desc.node, vl)) %>%
  mutate(vl1np1 = get_new_target(c("melUnd"), desc.node, vl)) %>%
  mutate(vl1np2 = get_new_target(c("calAnn"), desc.node, vl)) %>% as.data.table

dn.drop1.pval<- dn.drop %>% group_by(hog) %>% mutate(vl1p1p = compute_pval(dn.norm, vl1p1)) %>%
  mutate(vl1p2p = compute_pval(dn.norm, vl1p2)) %>%
  mutate(vl1p3p = compute_pval(dn.norm, vl1p3)) %>%
  mutate(vl1p4p = compute_pval(dn.norm, vl1p4)) %>%
  mutate(vl1p5p = compute_pval(dn.norm, vl1p5)) %>%
  mutate(vl1p6p = compute_pval(dn.norm, vl1p6)) %>%
  mutate(vl1np1p = compute_pval(dn.norm, vl1np1)) %>%
  mutate(vl1np2p = compute_pval(dn.norm, vl1np2)) %>%
  mutate(vlorigP = compute_pval(dn.norm, vl)) %>%
  select(hog, vl1p1p:vlorigP) %>% ungroup %>% distinct

dn.drop1.pval %>% select(-hog) %>% apply(., 2, function(x) qvalue(x)$pi0)

#need a programmatic way to construct next set of columns
vl_tips2 <- c("corBra", "ficAlb", "pseHum", "serCan", "taeGut", "geoFor", "melUnd", "calAnn", "none")
tipCount <- length(vl_tips2)
drop.res<-data.frame(dropped_tips = character(), q01 = numeric(), q05 = numeric(), p01 = numeric(), p001 = numeric(), pi0 = numeric())

for (i in 1:tipCount) {
  for (j in i:tipCount) {
    for (k in j:tipCount) {
    cols_to_drop = c(vl_tips2[i], vl_tips2[j], vl_tips2[k])
    colstring = paste0(cols_to_drop, collapse="-")
    ct<-dn.drop %>% group_by(hog) %>% mutate(newtarget = get_new_target(cols_to_drop, desc.node, vl), newpval = compute_pval(dn.norm, newtarget)) %>% ungroup %>% select(hog, newpval) %>% distinct %>% select(-hog) %>% apply(., 2, qvalue, pi0.method = "bootstrap")
    drop.res<-rbind(drop.res, data.frame(dropped_tips = colstring, q01 = sum(ct$newpval$qvalues <= 0.01), q05 = sum(ct$newpval$qvalues <= 0.05), p01 = sum(ct$newpval$pvalues <= 0.01), p001 = sum(ct$newpval$pvalues <= 0.001), pi0 = ct$newpval$pi0, row.names = NULL))
    print(drop.res)
   }
  }
}

#clean up
drop.res.clean <- drop.res %>% tbl_df %>% 
  mutate(calAnn = grepl("calAnn", dropped_tips), melUnd = grepl("melUnd", dropped_tips), corBra = grepl("corBra", dropped_tips), ficAlb = grepl("ficAlb", dropped_tips), pseHum = grepl("pseHum", dropped_tips), serCan = grepl("serCan", dropped_tips), taeGut = grepl("taeGut", dropped_tips), geoFor = grepl("geoFor", dropped_tips)) %>%
  mutate(os_drop_ct = as.numeric(corBra) + as.numeric(ficAlb) + as.numeric(pseHum) + as.numeric(serCan) + as.numeric(taeGut) + as.numeric(geoFor)) %>%
  mutate(nonos_drop_ct = as.numeric(calAnn) + as.numeric(melUnd)) %>% rowwise %>% 
  mutate(which_nonos = if(!calAnn & !melUnd) { "none" } else if (calAnn & melUnd) { "both" } else if (calAnn & !melUnd) { "calAnn" } else { "melUnd" }) %>%
  mutate(drop_ct = os_drop_ct + nonos_drop_ct) %>% select(-dropped_tips) %>% ungroup %>% distinct


drop.res.clean %>% ggplot(aes(as.factor(drop_ct),pi0)) + geom_jitter(aes(color=which_nonos), width = 0.1)

drop.res.clean %>% ggplot(aes(as.factor(drop_ct),q01)) + geom_jitter(aes(color=as.factor(os_drop_ct)), width = 0.1)

drop.res.clean %>% ggplot(aes(as.factor(drop_ct),p001)) + geom_jitter(aes(color=as.factor(which_nonos)), width = 0.15)

drop.res.clean %>% filter(drop_ct == 3) %>% mutate(non_osc = factor(which_nonos, levels=c('none', 'calAnn', 'melUnd', 'both'))) %>% lm(q01 ~ non_osc, data= . ) %>% summary

drop.res.clean %>% filter(drop_ct == 2) %>% mutate(non_osc = factor(which_nonos == "none")) %>% lm(q01 ~ non_osc, data= . ) %>% summary

write.table(drop.res.clean, file="drop_init_run.tsv", sep="\t", quote=F, row.names=F, col.names = T)

#FUNCTIONAL TESTS

dn.pval.merge <- dn.pval %>% left_join(., hog_to_gene)

#load Zhang data
zhang<-read.table("ZhangTableS28.txt", header=T, sep="\t", comment.char="", fill=T, stringsAsFactors = F) %>% tbl_df %>% select(gene = Gene.symbol, qval = adj.p.value)

dn.pval.merge <- dn.pval.merge %>% left_join(., zhang)
dn.pval.merge$vl.q = qvalue(dn.pval.merge$vl.p)$qvalue
dn.pval.merge <- dn.pval.merge %>% filter(!is.na(vl.q)) %>% mutate(zhang.sig05 = as.numeric(!is.na(qval)))

dn.pval.merge <- dn.pval.merge %>% left_join(., dn.dir)
dn.pval.merge$sig05 = as.numeric(dn.pval.merge$vl.q < 0.05) * dn.pval.merge$vl.dir

dn.pval.merge %>% with(., table(zhang.sig05, sig05 == 1)) %>% fisher.test
#significant if low magnitude overlap with Zhang et al 2014 vocal learning genes, when restricting to accelerated class (dir +)

#based on symbol matching so should be good but not perfect orthology


# PLOTTING BELOW ##

hog_to_plot = 11147
with(dn.default[dn.default$hog==hog_to_plot,], plot(sort(dn.norm), col=ifelse(ratite[order(dn.norm)],"red", "gray50"), pch=16))

#FIGURE 2A

#inset

plot(density(apply(dn.perm.pval, 2, function(x) sum(x < 0.05)/length(x))[2:length(dn.perm.pval)]), col="black", xlim=c(0,0.25), xlab="Fraction of tests with P < 0.05", las=1, bty="n", main="")
ratite.p.frac=sum(dn.pval$ratite.p[!is.na(dn.pval$ratite.p)] < 0.05)/sum(!is.na(dn.pval$ratite.p))
vl.p.frac=sum(dn.pval$vl.p[!is.na(dn.pval$vl.p)] < 0.05)/sum(!is.na(dn.pval$vl.p))
arrows(x0=vl.p.frac,y0=2,x1=vl.p.frac,y1=0, col="firebrick", lwd=2)
arrows(x0=ratite.p.frac,y0=2,x1=ratite.p.frac,y1=0, col="blue", lwd=2)
text(x=vl.p.frac, y=2.5, labels=c("Vocal Learners"))
text(x=ratite.p.frac, y=2.5, labels=c("Ratites"))

#figure -- need to do better than just overplotting 1000 lines though
permdist<-cut(as.data.frame(dn.perm.pval)[,2], breaks=seq(0,1,0.01),labels=F)
plot(table(permdist), type="l", col=rgb(100,100,100,alpha=50,maxColorValue=255), ylim=c(0,700), xaxt="n", ylab="Count", las=1, xlab="P-value")
for (i in 3:length(dn.perm.pval)) {
  permdist<-cut(as.data.frame(dn.perm.pval)[,i], breaks=seq(0,1,0.01))
  points(table(permdist), type="l", col=rgb(100,100,100,alpha=50,maxColorValue=255))
  
}

ratite<-cut(dn.pval$ratite.p, breaks=seq(0,1,0.01), labels=F)
vl<-cut(dn.pval$vl.p, breaks=seq(0,1,0.01), labels=F)

lines(table(ratite), type="l", col="blue", lwd=3, lty="dashed")
lines(table(vl), type="l", col="firebrick", lwd=3, lty="dashed")
axis(1, labels=seq(0,1,0.2), at=seq(0,100,20))
legend("topright", legend=c("random", "ratite", "vocal learners"), col=c("gray50", "blue", "firebrick"), lwd=3, lty="dashed")

