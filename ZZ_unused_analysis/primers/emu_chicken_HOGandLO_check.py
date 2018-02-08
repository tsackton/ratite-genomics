# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 17:39:24 2016

@author: Phil
"""

import csv
LOdic = {}
HOGdic = {}
fullList = []
LOlists = []
HOGlists = []
fout = open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/emuHOG_crossChecked.txt', 'w') #outfile
with open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/linkedCans.txt', 'rU') as handle: #infile
    reader=csv.reader(handle,delimiter='\t') 
    for strLine in reader: 
        fullList.append(strLine[2])
        fullSet=set(fullList)
with open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/emuMakers_perChickenHOG.txt', 'rU') as handle: #infile
    reader=csv.reader(handle,delimiter='\t') 
    for strLine in reader: 
        HOGlists.append([strLine[1],strLine[2]])
    for element in fullSet:
        for pair in HOGlists:
            if element in pair:
                HOGdic[element]=pair[0]
    for element in fullSet:
        if element not in HOGdic:
            HOGdic[element]="no HOG match"
with open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/droNov_allout_bed12out.txt', 'rU') as handle: #infile
    reader=csv.reader(handle,delimiter='\t') 
    for strLine in reader:
        LOlists.append([strLine[0],strLine[1]])        
    for element in fullSet:
        for pair in LOlists:
            if element in pair:
                LOdic[element]=pair[0]
    for element in fullSet:
        if element not in LOdic:
            LOdic[element]="no LO match" 
finalList=[]
for ID in LOdic:
    finalList.append([ID,LOdic[ID],HOGdic[ID]])
for ID in finalList:
    if ID[1]==ID[2]:
        ID.append('YES')
    else:
        ID.append('NO')

with open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/linkedCans.txt', 'rU') as handle: #infile
    reader=csv.reader(handle,delimiter='\t') 
    for strLine in reader:
        for ID in finalList:
            if strLine[2] in ID:
                if len(ID)<5:
                    ID.append(strLine[0])

for ID in finalList:
    fout.write('\t'.join(ID) +'\n')
fout.close()