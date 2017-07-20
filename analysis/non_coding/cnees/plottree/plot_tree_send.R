library(ape)
setwd("/Volumes/LaCie/Projects/Current/ratites/final/phyloAcc_May2017/")

## read in posterior of Z and loglik ##
postZ <- read.table("Combined_post_Z_06_2.txt", header=T)
score <- read.table("Combined_elem_lik2.txt", header=T)
score$BF2 <- -score$loglik_all + score$loglik_RES

setwd("~/Projects/birds/ratite_compgen/ratite-genomics/analysis/non_coding/cnees/")
#tree <- read.tree("hg38.phyloP100way.mam-a.mod")
tree <- read.tree("plottree/neut_ver3_final.named.mod")
species_readin <- read.table("plottree/species_names3.txt",stringsAsFactors = F)[,1]  ## order of species in the output files from our method

species <- tree$tip.label
S = length(species)

#target_species <- c("triMan1", "lepWed1", "odoRosDiv1", "orcOrc1", "turTru2") 
target_species <- c("strCam","rhePen","rheAme","casCas","droNov","aptRow","aptHaa","aptOwe","anoDid") 

tip.color = rep(1, length(species))
names(tip.color) <- species
tip.color[target_species] <- 4

## reorder the species same as the tree
node_idx <- 1:S
for(i in 1:(S-1)) # the last column edge from root ; the first 43 columns are species order the same in ape
{
  node_name <- tree$node.label[i]
  node_idx <-c(node_idx, which(species_readin == node_name))
}


## input: Z: one line from postZ, posterior of Z; tit: title of the figure (BF scores); offset=2, posterior of Z start from third column
## other input: node_idx, tree
## output: plot of one element
plotZPost <- function(Z,tit,offset=2) 
{
  
  ratio1 = Z[2] # c_rate
  ratio2 = Z[1] # n_rate
  
  tit = paste(tit, "\n rate1: ", round(ratio1,2), " rate2: ",round(ratio2,1))
  
  # color the tree by posterior mean of z
  E = length(Z)
  posterior_z <- rbind( Z[seq(offset + 2,E,by = 4)],Z[seq(offset +3,E,by = 4)], 
                        Z[seq(offset +4,E,by = 4)]) # get posterior of z 
  posterior_z <- posterior_z[,node_idx]
  
  
  posterior_z <- posterior_z[,tree$edge[,2]] # reorder z by edges
  
  mean_z <- colSums(posterior_z * 1:3);  
  
  rbPal <- colorRampPalette(c('gray','purple','green','red'))
  mytree= tree
  mytree$edge.length <-  tree$edge.length*(posterior_z[1,] + posterior_z[2,] * ratio1 + posterior_z[3,]* ratio2)
  
  edge_col <- rbPal(91)[round(mean_z*30)+1] #as.numeric(cut(posterior_z,breaks = 100))
  
  # color grey is missing
  missing = Z[seq(offset + 1,E,by = 4)]
  tip.color[missing >0] = "azure4"
  plot(mytree,edge.color = edge_col,tip.color = tip.color,edge.width  =2,label.offset = 0.01 ,main=tit)
  
}

mCE1623056
mCE225714
mCE300203
mCE1040658
mCE1084527
mCE1472811
mCE1095316
mCE1747138
mCE1753856
mCE967994 

## plot one element ##
id="mCE390119"
lk <- score %>% filter(ID==id)
tit = paste("BF1:", round(lk$log_ratio), "BF2:",round(lk$BF2), "Name:", lk$ID)
plotZPost(unlist(postZ[which(score$ID == id),]), tit, offset=2);

