# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 17:28:23 2016

@author: Phil
"""

import csv
ncbiList = []
hogList= []
fout = open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/emuMakers_perChickenHOG.txt', 'w') #outfile
with open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/linkedCans.txt', 'rU') as handle: #infile
    reader=csv.reader(handle,delimiter='\t') 
    for strLine in reader: 
        ncbiList.append(strLine[2])
with open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/new_hog_list.txt', 'rU') as handle: #infile
    reader=csv.reader(handle,delimiter='\t')
    hogLOL = []
    for strLine in reader:
        if strLine[3] == "galGal":
            if strLine[2] in ncbiList:
                hogLOL.append([strLine[0],strLine[2]]) #every pair in listoflists(LOL) is a list with two elements, a HOG and associated NCBI#
                hogList.append(strLine[0]) #this is a list of all chicken HOGs in ncbiList
with open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/new_hog_list.txt', 'rU') as handle: #infile
    reader=csv.reader(handle,delimiter='\t')
    for strLine in reader:
        if strLine[3] == "droNov":
            if strLine[0] in hogList:
                for pair in hogLOL:
                    if strLine[0] in pair:
                        fout.write(strLine[0] +'\t' + strLine[2] + '\t' + pair[1] + '\n')
fout.close()