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

from branchsite_parser import parse_hogs

def print_results(results, handle, model):
    for hog in results.keys():
        tree_len = len(results[hog]['trees'])
        res_len = len(results[hog]['results'])
        for i in range(1,res_len+1):
            trees = results[hog]['trees'].get(i, {})
            res = results[hog]['results'].get(i, {})
            #parse results
            for modelnum in sorted(res.get('NSsites',{}).keys(), key=int):
                #fixed options
                res_lnl = res.get('NSsites', {}).get(modelnum, {}).get('lnL')
                res_treelen = res.get('NSsites', {}).get(modelnum, {}).get('tree length')
                res_kappa = res.get('NSsites', {}).get(modelnum, {}).get('parameters', {}).get('kappa')
                #what to get depends on the model number
                if modelnum == 0:
                    #m0 case
                    res_omega = res.get('NSsites', {}).get(modelnum, {}).get('parameters', {}).get('omega')
                elif modelnum == 1 or modelnum == 2:
                    #m1/m2 case, use site classes 
                    res_omega = None
                    res_omegas = res.get('NSsites', {}).get(modelnum, {}).get('parameters', {}).get('site classes', {})
                    for sc in sorted(res_omegas.keys(), key=int):
                        if res_omega is None:
                            res_omega = str(round(res_omegas[sc]['proportion'],3)) + ":" + res_omegas[sc]['omega']
                        else:
                            res_omega += "," + str(round(res_omegas[sc]['proportion'],3)) + ":" + res_omegas[sc]['omega']
                elif modelnum == 7:
                    #m7 case
                    m7_p = str(round(float(res.get('NSsites', {}).get(modelnum, {}).get('parameters', {}).get('p')),3))
                    m7_q = str(round(float(res.get('NSsites', {}).get(modelnum, {}).get('parameters', {}).get('q')),3))
                    res_omega = m7_p + "/" + m7_q
                elif modelnum == 8:
                    m8_p = str(round(float(res.get('NSsites', {}).get(modelnum, {}).get('parameters', {}).get('p')),3))
                    m8_q = str(round(float(res.get('NSsites', {}).get(modelnum, {}).get('parameters', {}).get('q')),3))
                    m8_p0 = str(round(float(res.get('NSsites', {}).get(modelnum, {}).get('parameters', {}).get('p0')),3))
                    m8_p1 = str(round(1 - float(m8_p0),3))
                    m8_w = str(round(float(res.get('NSsites', {}).get(modelnum, {}).get('parameters', {}).get('w')),3))
                    res_omega = m8_p0 + ":" + m8_p + "/" + m8_q + "," + m8_p1 + ":" + m8_w
                else:
                    pass
                    
                print(hog, model, i, trees.get('is_species_tree'), trees.get('original'), modelnum, res_lnl, res_treelen, res_kappa, res_omega, sep="\t", end="\n", file=handle)
 

def main():
    print("Starting to parse.")
    hogfile_toparse = sys.argv[1]
    resfile_foroutput = sys.argv[2]
    with open(hogfile_toparse) as hfile:
        hogs=[line.rstrip("\n") for line in hfile]
    
    print("Done getting files.")
    with open(resfile_foroutput, 'w') as ofile:
        print("hog", "model", "treenum", "species_tree", "newick_string", "model_num", "lnl", "treelen", "kappa", "omega", sep="\t", end="\n", file=ofile, flush=True)
        results = parse_hogs(hogs,"site",["paml"],True,True)
        print_results(results, ofile, "site")

if __name__ == "__main__":
    main()
