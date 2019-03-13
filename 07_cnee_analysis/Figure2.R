######## draw acceleration tree and sequence alignment ########
library(seqinr)
library(ggplot2)
library(reshape2)
library(ape)


#### load phylogenetic tree and common name of species ####
## input: tree_path is file of phylogenetic tree; 
#        species_name is a file output by PhyloAcc containing the species name for each column in *postZ* files; 
#        common_name is optional file with three columns: abbreviation of species name appeared in data files and output files; full species names; and species comman name shown on the plot 
## output: a list as input to plotZPost and plotAlign to generate plots
prepare_data <- function(tree_path = "neut_original.named.mod", species_name = "species_names.txt", common_name = "birdname2.txt")
{
  
  tree <- read.tree(tree_path) # 43 species
  species <- tree$tip.label
  
  ## reorder species in the output file to be aligned with tree
  species_readin <- read.table(species_name,stringsAsFactors = F)[,1]
  
  node_idx <- c()
  for(node_name in tree$tip.label) # the last column edge from root ; the first 43 columns are species order the same in ape
  {
    node_idx <-c(node_idx, which(species_readin == node_name))
  }
  
  for(node_name in tree$node.label) # the last column edge from root ; the first 43 columns are species order the same in ape
  {
    node_idx <-c(node_idx, which(species_readin == node_name))
  }
  
  ## species common names
  if(!is.null(common_name))
  {
    comnam <- read.table(common_name,sep="\t", row.names = 1, stringsAsFactors = F)
    return(list("tree" = tree, "node_idx" = node_idx, "tip" = comnam[tree$tip.label,2]))
  }
  
  return(list("tree" = tree, "node_idx" = node_idx, "tip" = tree$tip.label))
  
}

#### Plot the acceleration pattern (posterior of substitution rate on each branch) of one element ####
## input: Z: one row from *postZ* file, posterior of Z for each branch; 0: missing(only for outgroup), 1: neutral, 2: conserved, 3: accelerated
#         tit: title of the figure (e.g. BF scores); 
#         treeData: output from  prepare_data
#         target_species: group of species shown in different color (e.g. phenotypically convergent species), please use the abbrev. name same as the tree
#         offset=2, posterior of Z start from the third column
## output: plot of one element

plotZPost <- function(Z, treeData, target_species=NULL, tit=NULL, idx_offset=3, cex.score = 2)
{
  
  species <- treeData$tree$tip.label
  tip.color = rep(1, length(species))
  names(tip.color) <- species
  if(!is.null(target_species)) tip.color[target_species] <- 4
  
  ratio1 = Z[2] # c_rate #offset-3
  ratio2 = Z[1] # n_rate #offset-4
  #print(ratio1);print(ratio2);
  
  
  # color the tree by posterior mean of z
  E = length(Z)
  posterior_z <- rbind( Z[seq(idx_offset+2,E,by = 4)],Z[seq(idx_offset+3,E,by = 4)], 
                        Z[seq(idx_offset+4,E,by = 4)]) # get posterior of z (only 1???2???3)
  
  # get number of losses by posterior of Z
  #loss = sum(posterior_z[3, 1:length(species)]) - sum(posterior_z[3, -1:-length(species)])
  
  posterior_z <- posterior_z[,treeData$node_idx]
  posterior_z <- posterior_z[,treeData$tree$edge[,2]] # reorder z by edges
  
  mean_z <- colSums(posterior_z * 1:3);  
  
  rbPal <- colorRampPalette(c('gray','goldenrod','navyblue','firebrick1'))
  mytree= treeData$tree
  mytree$edge.length <-  mytree$edge.length*(posterior_z[1,] + posterior_z[2,] * ratio1 + posterior_z[3,]* ratio2)
  mytree$tip.label <- treeData$tip
  edge_col <- rbPal(91)[round(mean_z*30)+1] 
  
  # color grey is missing
  missing = Z[seq(idx_offset + 1,E,by = 4)]
  tip.color[missing >0] = "azure4"
  par(mar= c(2,0,0,0))
  plot(mytree,edge.color = edge_col,tip.color = tip.color,edge.width  =3,label.offset = 0.003,no.margin =F, cex=1.1, bg=NA)
  mtext(substitute(paste(t, r[1], "=", r1, ", ", r[2], "=", r2),
                   list(t=tit, r1 =round(ratio1,2), r2 = round(ratio2,2))), side = 1, cex=cex.score, line = 0.5)
  
}

## read in scores and posterior of Z ##
data_path <- "/Users/tim/Projects/birds/ratite_compgen/DRYAD/07_cnees/results/orig_v2_phyloAcc"
version <- "gain"
score <- read.table(paste0(data_path, "/", version, "/", "Combined_elem_lik.txt"), header=T)
treeData <- prepare_data(tree_path = "/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/neut_original.named.mod",
                         species_name = "/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/species_names.txt",
                         common_name = "/Users/tim/Projects/birds/ratite_compgen/ratite-genomics/07_cnee_analysis/birdname2.txt")

treeData$target_species <- c("strCam","rhePen","rheAme","casCas","droNov","aptRow","aptHaa","aptOwe","anoDid") 

postZ <- read.table(paste0(data_path, "/", version, "/", "Combined_post_Z_M2.txt"), header=T)


# selected element to plot
elements <- c("mCE1333824", "mCE967994", "mCE1140641") 
ks = which(score$ID %in% elements)-1

save_to_pdf <- TRUE

for(k in ks) 
{
  if (save_to_pdf) {
    pdf(paste("tree_",version,"_",score$ID[k+1],".pdf",sep=""))
  }
  lk <- subset(score,No.==k); print(lk); print(postZ[k+1,1:2])
  tit = paste("logBF1:", round(lk$logBF1), "logBF2:",round(lk$logBF2), "  ")
  plotZPost(unlist(postZ[k+1,]),target_species=treeData$target_species, tit=tit, treeData,idx_offset=5);

  if (save_to_pdf) {
    dev.off()
  }
}

