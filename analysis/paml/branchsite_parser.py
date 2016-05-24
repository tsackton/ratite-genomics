#!/usr/bin/env python

import os
import sys
from Bio import Phylo
from Bio.Phylo.PAML import codeml
from Bio.Phylo.PAML import _parse_codeml
import io
import re

def _bitstrs(tree):
    #function from http://biopython.org/wiki/Phylo
    # store and return all _BitStrings    
    bitstrs = set()
    term_names = [term.name for term in tree.get_terminals()]
    term_names.sort()
    for clade in tree.get_nonterminals():
        clade_term_names = [term.name for term in clade.get_terminals()]
        boolvals = [name in clade_term_names for name in term_names]
        bitstr = _BitString(''.join(map(str, map(int, boolvals))))
        bitstrs.add(bitstr)
    return bitstrs

def compare_trees(tree1, tree2):
    #function from http://biopython.org/wiki/Phylo
    term_names1 = [term.name for term in tree1.get_terminals()]
    term_names2 = [term.name for term in tree2.get_terminals()]
    # false if terminals are not the same
    if set(term_names1) != set(term_names2):
        return False
    # true if _BitStrings are the same
    if _bitstrs(tree1) == _bitstrs(tree2):
        return True
    else:
        return False

def classify_tree (tree):
    #take a tree object, identify the clade(s) that are the foreground branches
    #return tree object with #1 removed and the foreground branches IDed
    
    regex = re.compile(r"#")
    foreground = []
    
    for term in tree.get_terminals():
        if regex.search(term.name):
            term.name = term.name.replace("#1","")
            foreground.append(term.name)
    
    return (tree, foreground)

def get_trees (file,speciestree):
    #read a tree file, strip out the spaces from the individual trees, and read each one
    #returns list of trees where list index + 1 is the tree number
    
    with open(file, 'r') as tf:
         tree_strings = tf.readlines()
    
    paml_trees = {}
    #tree strings need to get converted to trees now
    for i in range(1,len(tree_strings)):
         cur_tree = tree_strings[i].replace(" ","")
         processed_trees = Phylo.read(io.StringIO(cur_tree), "newick")
         final_tree, foreground_sp = classify_tree(processed_tree)
         sptree = compare_trees(final_tree, speciestree)
         tree_dict = {'tree':final_tree, 'foreground':foreground_sp, 'is_species_tree':sptree}
         paml_trees[i] = tree_dict
        
    return(paml_trees)

def parse_codeml_string (handle):
    results = {}
    lines = handle.readlines()
    (results, multi_models, multi_genes) = _parse_codeml.parse_basics(lines, results)
    results = _parse_codeml.parse_nssites(lines, results, multi_models, multi_genes) 
    results = _parse_codeml.parse_pairwise(lines, results)
    results = _parse_codeml.parse_distances(lines, results)
    if len(results) == 0:
        raise ValueError("Invalid results file") 
    return results 
     
def split_results (file):
     #takes a paml output file, splits on TREE #, returns list of output files to parse
     
    with open(file, 'r') as pf:
        results=pf.read()
    
    #split results
    pieces = results.split("\nTREE")
    header = pieces.pop(0)
    paml_results = {}
    for i in range(len(pieces)):
        recons_res = header + "\nTREE" + pieces[i]
        parsed = parse_codeml_string(io.StringIO(recons_res))
        tree_index = i+1
        paml_results[tree_index] = parsed
    return(paml_results)

if __name__ == "__main__":
    
    