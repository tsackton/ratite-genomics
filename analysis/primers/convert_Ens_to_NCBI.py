# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 13:35:26 2016

@author: Phil
"""
import csv

class Vividict(dict): #used for autovivification of dictionary 
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value
        
cdic = Vividict() #make dictionary
fout = open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/linkedCans.txt', 'w') #outfile
with open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/ce_annotation.tsv', 'rU') as handle: #infile
    reader=csv.reader(handle,delimiter='\t') 
    for strLine in reader:
        cdic[strLine[7]][strLine[0]]=strLine[9] #7 is ensembl ID, 0 is mCE number, 9 is NCBI accession     
with open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/ceCans.txt', 'rU') as handle: #infile
    reader=csv.reader(handle,delimiter='\t') 
    for strLine in reader: #using mCE numbers from first candidate set...
        for EnsemblID in cdic.keys(): #for each ensembl ID
            for mCE in cdic[EnsemblID].keys(): #for each mCE and NCBI associated with the ensembl ID
                if strLine[0] == mCE: #if our candidate mCE is present
                    fout.write(EnsemblID + '\t' + strLine[0] + '\t' + cdic[EnsemblID][mCE] + '\n')                    
fout.close()