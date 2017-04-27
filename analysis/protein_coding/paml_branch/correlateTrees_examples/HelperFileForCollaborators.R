# Load required modules
source("ERCfunctions161214.R")


# Create a "trees object" containing a single tree for each gene.
#	Input file format is "name, tab, Newick tree string, newline".
#	Species names must be consistent throughout.
trees=readTrees("20DRO_aaML.FILTERED.txt", rearrange=T, computePaths=T)
save(trees, file="droso_tree_obj.RData")	# Save as R data for rapid loading next time.


### Correlate trees with binary traits ###
# read in the tree representing the binary trait to be correlated. Branches with 1 are foreground, 0 is background.
trait_tree <- read.tree("bin_marine.tre")

# compute the correlation b/w trait and each tree
result <- correlateTreesBinary(trees , trait_tree)

	# The result is a list of 2 vectors, "r" and "p". The values are in the same order as the rownames of the trees$trees object.
	# "r" contains the AUC metric for each region.
		# r can be thought of as the proportion of background branches for which the foreground branches have a higher evolutionary rate.
		# r > 0.5 means the foreground is more rapid than the background on average
		# r < 0.5 means the background is more rapid
	# "p" contains the p-values for these same relationships




### Correlate trees with continuous traits ###
# THIS FUNCTION IS BEING COMPLETELY REVAMPED TO BE MORE ROBUST TO VARIABLE BRANCH LENGTHS
# NEW FUNCTION COMING SOON




# Compute full ERC matrix
#	Testing may be done with the "maxDo" argument to limit the number of values to compute.
#	species.threshold sets a minimum number of species to be SHARED between a pair of genes, otherwise no ERC value is computed.
#	species.list allows user to study a subset of species
load("droso_tree_obj.RData")
#erc=correlateTreesAll(trees,  usePaths=T, maxDo=10)   # Stops after 10 ERC values are calculated for troubleshooting
erc=correlateTreesAll(trees,  usePaths=T, species.threshold=10, species.list=c("DGRI","DVIR","DMOJ","DBIP","DANA","DKIK","DFIC","DMEL","DSEC","DYAK","DERE","DEUG","DBIA","DTAK","DELE","DRHO","DPSE","DWIL"))
save(erc, file="FILENAME");		# Save ERC matrix for R
write.table( t(erc), "erc_droso20.tsv", sep="\t", row.names=F, col.names=T , quote=F)		# writes tab-delimited table for other applications

# small function to clean out empty rows and columns from the ERC matrix resulting when genes have an insufficient number of species.
cerc=removeEmpties(erc)


# Code to study a correlation between a specified pair of genes
# This will give you enough data to compute the ERC value for a specified pair instead of the full genome-wide matrix.
#	The "res" is a projection object containing the original branch values "l1 l2", the normalized branch-specific rates "e1 e2", the normalization vector "nv", etc...
res=correlateTreesProj("gene1", "gene2", trees, plot=T, usePaths = T)
#	This plots the normalized values of gene1 versus gene2 to see the rates that led to the ERC value.
plotWtext(res$e1, res$e2, labels=res$names, cex=0.7)
