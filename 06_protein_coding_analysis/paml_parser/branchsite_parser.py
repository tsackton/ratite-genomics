#!/usr/bin/env python

import os
import sys
from Bio import Phylo
from Bio.Phylo.PAML import codeml
from Bio.Phylo.PAML import _parse_codeml
from Bio.Phylo.Consensus import _BitString
import io
import re

def _bitstrs(tree,term_names):
    #function from http://biopython.org/wiki/Phylo
    # store and return all _BitStrings    
    bitstrs = set()
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
    if _bitstrs(tree1,term_names1) == _bitstrs(tree2,term_names2):
        return True
    else:
        return False

def classify_tree (tree):
    #take a tree object, identify the clade(s) that are the foreground branches
    #return tree object with #1 removed and the foreground branches IDed
    
    regex = re.compile(r"#")
    foreground = []
    genetree_names = re.compile(r"([A-Za-z])_[A-Za-z0-9_.-]+")    
    for term in tree.get_terminals():
        if regex.search(term.name):
            term.name = term.name.replace("#1","")
            if len(term.name) > 6:
                genetree_match = genetree_names.search(term.name)
                if genetree_match:
                    term.name = genetree_match.group(1)
            foreground.append(term.name)
    
    foreground_string = ":".join(sorted(foreground))
    return (tree, foreground_string)

def parse_trees (file,speciestree):
    #read a tree file, strip out the spaces from the individual trees, and read each one
    #returns list of trees where list index + 1 is the tree number
    
    with open(file, 'r') as tf:
         tree_strings = tf.readlines()
    
    paml_trees = {}
    seen_trees = {}
    #tree strings need to get converted to trees now
    for i in range(1,len(tree_strings)):
        if tree_strings[i] in seen_trees:
            treekey = seen_trees[tree_strings[i]]
            paml_trees[i] = paml_trees[treekey]
        else:
            seen_trees[tree_strings[i]] = i
            cur_tree = tree_strings[i].replace(" ","")
            processed_trees = Phylo.read(io.StringIO(cur_tree), "newick")
            final_tree, foreground_sp = classify_tree(processed_trees)
            if speciestree is not None:
                sptree = compare_trees(final_tree, speciestree)
            else:
                sptree = False
            tree_dict = {'tree':final_tree, 'foreground':foreground_sp, 'is_species_tree':sptree, 'original':tree_strings[i].rstrip()}
            paml_trees[i] = tree_dict
    
    return(paml_trees)

def parse_codeml_string (handle):
    results = {}
    lines = handle.readlines()
    (results, multi_models, multi_genes) = _parse_codeml.parse_basics(lines, results)
    results = _parse_codeml.parse_nssites(lines, results, multi_models, multi_genes) 
#    results = _parse_codeml.parse_pairwise(lines, results)
#    results = _parse_codeml.parse_distances(lines, results)
    if len(results) == 0:
        raise ValueError("Invalid results file") 
    return results 
     
def parse_multitree_results (file):
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

def parse_multitree_multimodel_results (file):
     #takes a paml output file, splits on TREE #, returns list of output files to parse
     
    with open(file, 'r') as pf:
        results=pf.read()
    
    #split results
    models = results.split("\nModel ")
    header = models.pop(0)
    paml_strings = {}
    paml_results = {}
    for model in models:
        pieces = model.split("\nTREE")
        subheader = "\nModel " + pieces.pop(0)
        for i in range(len(pieces)):
            tree_index = i+1
            try:
                paml_strings[tree_index] += subheader + "\nTREE" + pieces[i]
            except KeyError:
                paml_strings[tree_index] = header + subheader + "\nTREE" + pieces[i]
    for tree_num in paml_strings:
        paml_results[tree_num] = parse_codeml_string(io.StringIO(paml_strings[tree_num]))

    return(paml_results)

def parse_hogs(hoglist,model,basedir,verbose=True,multisite=False):
    #take list of hogs, return parsed final results dictionary
    final_results = {}
    for hog in hoglist:
        if verbose:
            print("Working on", hog, flush=True)
            
        toppath = '{:0>4}'.format(int(hog) % 100)
        # 0000/100/100.codeml.ancrec.ctl.out/
        fullpath = basedir + "/" + toppath + "/" + hog + "/" + hog + ".codeml." + model + ".ctl.out"
        results_file = fullpath + "/" + model + ".out"
        control_file = fullpath + "/" + hog + ".codeml." + model + ".ctl"
        #get species tree
        sptreepath = basedir + "/" + toppath + "/" + hog + "/" + hog + ".final_spt.nwk"
        try:
            species_tree = Phylo.read(sptreepath, "newick")
        except FileNotFoundError:
            species_tree = None
                
        cml = codeml.Codeml()
        try:
            cml.read_ctl_file(control_file)
        except OSError:
            print("Couldn't parse file for", hog, "at", pamldir + "/" + fullpath)
            continue
            
        tree_file = fullpath + "/" + cml.tree
        #now process
        parsed_trees = parse_trees(tree_file,species_tree)
        try:
            if multisite:
                parsed_results = parse_multitree_multimodel_results(results_file)
            else:
                parsed_results = parse_multitree_results(results_file)
        except FileNotFoundError:
            print("Couldn't parse file for", hog, "at", pamldir + "/" + fullpath)
            continue
            
        #check that we have a result for each tree
        if len(parsed_trees) < len(parsed_results):
            print("Warning, too few trees for number of results for", hog, "in", results_file)
            continue
        elif len(parsed_trees) > len(parsed_results):
            #remove trees that aren't in results
            trimmed_trees = {x:parsed_trees[x] for x in parsed_results.keys()}
            parsed_trees = trimmed_trees
                
        if hog in final_results:
            #append
            cur_len = len(final_results[hog]['trees'])
            if cur_len != len(final_results[hog]['results']):
                print("Warning, something went wrong!!")
                    
            #update keys (tree numbers)
            new_trees = {int(x)+cur_len:parsed_trees[x] for x in parsed_trees.keys()}
            new_results = {int(x)+cur_len:parsed_results[x] for x in parsed_results.keys()}
            final_results[hog]['trees'].update(new_trees)
            final_results[hog]['results'].update(new_results)
                    
        else:       
            final_results[hog] = {'trees':parsed_trees, 'results':parsed_results}
    
    return(final_results) 
    

def print_results (results, handle, model):
    #results is a complicated dict, but the basic format is: hog, tree->tree_dict, results->results_dict
    #tree_dict and results_dict are matched by key
    #to test let's start by printing out: hog id, tree num, foreground, species_tree vs gene_tree
    #original tree, then paml parameters
    for hog in results.keys():
        tree_len = len(results[hog]['trees'])
        res_len = len(results[hog]['results'])
#       print(tree_len, res_len)
#       print(results[hog]['trees'].keys())
#       print(results[hog]['results'].keys())
        for i in range(1,res_len+1):
            trees = results[hog]['trees'].get(i, {})
            res = results[hog]['results'].get(i, {})
            print(hog, model, i, trees.get('foreground'), trees.get('is_species_tree'), trees.get('original'), sep="\t", end="\t", file=handle)
            #parse results
            res_lnl = res.get('NSsites', {}).get(2, {}).get('lnL')
            res_siteclass = res.get('NSsites', {}).get(2, {}).get('parameters', {}).get('site classes')
            res_treelen = res.get('NSsites', {}).get(2, {}).get('tree length')
            print(res_lnl, res_treelen, sep="\t", end="", file=handle)
            #need to iterate over site class numbers
            try:
                for sc in sorted(res_siteclass):
                    print("\t" + res_siteclass[sc]['proportion'], res_siteclass[sc]['branch types']['foreground'], res_siteclass[sc]['branch types']['background'], sep="\t", end = "", file=handle)
                print(file=handle)
            except:
                print("\tNone", 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None', sep="\t", end="\n", file=handle)

def main():

    print("Starting to parse.")
    hogfile_toparse = sys.argv[1]   
    resfile_foroutput = sys.argv[2] 
    #get hog list, for now from a file
    with open(hogfile_toparse) as hfile:
        hogs=[line.rstrip('\n') for line in hfile]

    print("Done getting files.", flush=True)                
    with open(resfile_foroutput, 'w') as ofile:
        print("hog", "model", "treenum", "foreground_species", "species_tree", "newick_string", "lnl", "treelen", "class0_prop", "class0_fore", "class0_back", "class1_prop", "class1_fore", "class1_back", "class2a_prop", "class2a_fore", "class2a_back", "class2b_prop", "class2b_fore", "class2b_back", sep="\t", end="\n", file=ofile, flush=True)

        for model in ("branchsite", "branchsitenull"):
            results = parse_hogs(hogs,model,["paml", "paml_branch"])
            #print out
            print_results(results, ofile, model)
                 
if __name__ == "__main__":
    main()
          
            
            
            