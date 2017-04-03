# -*- coding: utf-8 -*-
"""
Created on Sat Apr  1 18:14:45 2017

@author: Phil
"""
import csv

fileout=('/Users/Phil/Dropbox/Doctorate/atac/MarchPeakAnalysis/woForNewBayesCategories_DFnames.txt')
fout= open(fileout, 'w') 

DFnames = ["FLtxtlooseBioGeo","FLtxtlooseDollo","FLtxtstrictBioGeo","FLtxtstrictDollo","HLtxtlooseBioGeo","HLtxtlooseDollo","HLtxtstrictBioGeo","HLtxtstrictDollo","Keel10txtlooseBioGeo","Keel10txtlooseDollo","Keel10txtstrictBioGeo","Keel10txtstrictDollo","Keel9txtlooseBioGeo","Keel9txtlooseDollo","Keel9txtstrictBioGeo","Keel9txtstrictDollo","Pec10txtlooseBioGeo","Pec10txtlooseDollo","Pec10txtstrictBioGeo","Pec10txtstrictDollo","Pec9txtlooseBioGeo","Pec9txtlooseDollo","Pec9txtstrictBioGeo","Pec9txtstrictDollo","Strn10txtlooseBioGeo","Strn10txtlooseDollo","Strn10txtstrictBioGeo","Strn10txtstrictDollo","Strn9txtlooseBioGeo","Strn9txtlooseDollo","Strn9txtstrictBioGeo","Strn9txtstrictDollo"]
DFmods = []
for element in DFnames:
    DFmods.append(element.split("txt"))
    
with open('/Users/Phil/Dropbox/Doctorate/atac/MarchPeakAnalysis/woForNewBayesCategories.txt', 'rb') as handle:
    reader=csv.reader(handle,delimiter=' ')
    for strLine in reader:
        for name in DFnames:
            if name.split("txt")[0] in strLine[0] and name.split("txt")[1] in strLine[0]:
                strLine[0] = name
                fout.write("\t".join(strLine)+"\n")
fout.close()