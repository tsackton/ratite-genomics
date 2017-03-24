#!/usr/bin/env python

import os
import sys
from Bio import Phylo
from Bio.Phylo.PAML import codeml
from Bio.Phylo.PAML import _parse_codeml
from Bio.Phylo.Consensus import _BitString
import io
import re

#should refactor this into an actual module I guess

from branchsite_parser import _bitstrs
from branchsite_parser import compare_trees
from branchsite_parser import classify_tree
from branchsite_parser import parse_trees
from branchsite_parser import parse_codeml_string
from branchsite_parser import parse_multitree_results
from branchsite_parser import parse_hogs

def parse_rst_mutations (file):
    #takes a path to an rst file and returns a parsed version as a dict
    #going to be complicated....   
    #print("Starting to parse rst file.")
    mut_results = {}
    #define some starting variables
    cur_tree = 0
    cur_branch = None
    #compile regexes
    treematch = re.compile('^TREE\s*\#\s*(\d+)') #group 1 is tree number
    branchmatch = re.compile('^Branch\s*(\d+):\s*(\d+\.\.\d+)') #group 1 is branch number, group 2 is branch key
    mutmatch = re.compile('^\s+(\d+)\s+\w+\s+\((\w)\).+\((\w)\)') #group 1 is position, group 2 is start mut, group 3 is end mut
    with open(file, 'r') as rst:
        for line in rst:
            treetest = treematch.match(line)
            branchtest = branchmatch.match(line)
            muttest = mutmatch.match(line)            
            if treetest:
                cur_tree = treetest.group(1)
                #print("Got match at ", line, flush=True)
                #print("Tree is now ", cur_tree, flush=True)
                mut_results[cur_tree] = {}            
            if branchtest:
                #print("Got match at ", line, flush=True)
                cur_branch = branchtest.group(2)           
                #print("Branch is now ", cur_branch, flush=True)
            if muttest:
                pos = muttest.group(1)
                fromaa = muttest.group(2)
                toaa = muttest.group(3)
                if fromaa == toaa:
                     continue
                event = fromaa + "-" + toaa       
                if pos not in mut_results[cur_tree]:
                    mut_results[cur_tree][pos] = {}
                
                mut_results[cur_tree][pos][cur_branch] = event         
    return(mut_results)

def print_results(results, handle, model):
    for hog in results.keys():
        tree_len = len(results[hog]['trees'])
        res_len = len(results[hog]['results'])
        for i in range(1,res_len+1):
            trees = results[hog]['trees'].get(i, {})
            res = results[hog]['results'].get(i, {})
            print(hog, model, i, trees.get('is_species_tree'), trees.get('original'), sep="\t", end="\t", file=handle)
            #parse results
            res_lnl = res.get('NSsites', {}).get(0, {}).get('lnL')
            res_omega = res.get('NSsites', {}).get(0, {}).get('parameters', {}).get('omega')
            res_treelen = res.get('NSsites', {}).get(0, {}).get('tree length')
            res_kappa = res.get('NSsites', {}).get(0, {}).get('parameters', {}).get('kappa')
            print(res_lnl, res_treelen, res_kappa, res_omega, sep="\t", end="\n", file=handle)

def parse_hogs_for_rst(hoglist,model,basedir,verbose=True):
    #take list of hogs, return parsed final results dictionary
    final_results = {}
    for hog in hoglist:
        if verbose:
            print("Working on ", hog, flush=True)
            
        toppath = '{:0>4}'.format(int(hog) % 100)
        # 0000/100/100.codeml.ancrec.ctl.out/
        fullpath = basedir + "/" + toppath + "/" + hog + "/" + hog + ".codeml." + model + ".ctl.out"
        results_file = fullpath + "/" + "rst"
        try:
            single_results = parse_rst_mutations(results_file)
            final_results[hog] = single_results
        except FileNotFoundError:
            continue
                
    return(final_results)

def print_rst_results(results,handle):
    for hog in results.keys():
        for tree in results[hog].keys():
            for pos in sorted(results[hog][tree].keys(), key=int):
                events = ','.join('{}:{}'.format(k,v) for (k,v) in results[hog][tree][pos].items())
                print(hog, tree, pos, events, sep="\t", end="\n", file=handle)
            

def main():
    print("Starting to parse.")
    hogfile_toparse = sys.argv[1]
    resfile_foroutput = sys.argv[2]
    mutfile_foroutput = sys.argv[3]
    basedir = sys.argv[4]
    with open(hogfile_toparse) as hfile:
        hogs=[line.rstrip("\n") for line in hfile]
    
    print("Done getting files.")
    with open(resfile_foroutput, 'w') as ofile:
        print("hog", "model", "treenum", "species_tree", "newick_string", "lnl", "treelen", "kappa", "omega", sep="\t", end="\n", file=ofile, flush=True)
        results = parse_hogs(hogs,"ancrec",basedir)        
        print_results(results, ofile, "ancrec")
    
    print("Done parsing model output, starting on ancestral reconstructions.")
    with open(mutfile_foroutput, 'w') as mutfile:
        print("hog", "tree", "pos", "mutations", sep="\t", end="\n", file=mutfile, flush=True)
        rst_results = parse_hogs_for_rst(hogs,"ancrec",basedir)       
        print_rst_results(rst_results, mutfile)

if __name__ == "__main__":
    main()
