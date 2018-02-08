#! /usr/bin/env python

import sys
import os
import Bio
from Bio import Phylo
from io import StringIO
from StringIO import StringIO
import io
import re


#2 arguments required = PAML output followed by HOG ID.
try:
	paml_output_file = sys.argv[1]
	tree_1_output_file = sys.argv[2]+".tree1.nwk"
	tree_2_output_file = sys.argv[2]+".tree2.nwk"
except:
	print "PAML output file name and HOG ID required as input"

paml_output = open(paml_output_file,"r")
paml_results = paml_output.read()



#Extract all trees from a list of text that have a six letter abbreviation following the pattern: genSpe. Returns a list of trees.
def extract_trees(list):
	searchRegex = re.compile("[a-z]{3}[A-Z]{1}[a-z]{2}").search
	trees = [l for l in list for m in [searchRegex(l)] if m]
	return trees
	
#Divide all branch lengths of newick tree by 1/3
def change_brlens(tree):
	tree_sp1 = tree.split(":")
	newtree = []
	regex = re.compile("[0-9]{1}\.[0-9]{1,6}")
	for group in tree_sp1:
		if regex.search(group):
			#print group
			if re.match("[0-9]{1}\.[0-9]{1,6}\)",group):
				group_spl = group.split(")")
				bl = float(group_spl[0])/3
				group_spl[0] = str(bl)
				bl = ")".join(group_spl)
				newtree.append(bl)
			elif re.match("[0-9]{1}\.[0-9]{1,6}\,",group):
				group_spl = group.split(",")
				bl = float(group_spl[0])/3
				group_spl[0] = str(bl)
				bl = ",".join(group_spl)
				newtree.append(bl)
			else:
				newtree.append(group)
		else:
			newtree.append(group)
		
	newtree = ":".join(newtree)
	return newtree


#Split results on "TREE", get TREE #  1 and TREE #  2. Look for (speSpe to ID string with tree newick and species names.
pieces = paml_results.split("\nTREE")
for piece in pieces:
	if piece[1:5] == "#  1":
		tree_1_output = open(tree_1_output_file,"w")
		tree1_all = piece.split("\n")
		tree1 = extract_trees(tree1_all)
		tree1 = tree1[0].replace(" ","")
		tree1_newBL = change_brlens(tree1)
		tree_1_output.write(tree1_newBL)
		tree1_newBL = Phylo.read(StringIO(tree1_newBL),"newick")
		#Phylo.draw_ascii(tree1_newBL)
		#tree_1_output.close()
		
	elif piece[1:5] == "#  2":
		tree_2_output = open(tree_2_output_file,"w")	
		tree2_all = piece.split("\n")
		tree2 = extract_trees(tree2_all)
		tree2 = tree2[0].replace(" ","")
		tree2_newBL = change_brlens(tree2)
		tree_2_output.write(tree2_newBL)
		#tree2_newBL = Phylo.read(StringIO(tree2_newBL),"newick")	
		#Phylo.draw_ascii(tree2_newBL)
		tree_2_output.close()

	else:
		pass
		

paml_output.close()
