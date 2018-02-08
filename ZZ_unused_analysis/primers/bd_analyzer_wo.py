# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 12:24:50 2015

@author: phil-grayson
"""

import re
import sys
import csv

phile = str(sys.argv[1])
fout = open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/'+phile.split('.')[0]+'out.txt', 'w') #outfile
with open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/'+phile, 'rU') as handle: #infile
    reader=csv.reader(handle,delimiter='\t') 
    mainlist = []    
    for strLine in reader: 
        metaData1 = strLine[3] #metaData contains info from chicken - gene ID, genebank ID, start and stop of exon
        GeneIDsearch = re.search(r'GeneID\:([0-9]*)\,', metaData1)           
        geneID = GeneIDsearch.group(1)        
        metaData2 = strLine[15]
        makerIDsearch = re.search(r'ID\=(DNOV[0-9]*)\-', metaData2)           
        makerID = makerIDsearch.group(1)
        mainlist.append([makerID,geneID])
    mainset = set(map(tuple,mainlist))
    newlist = map(list,mainset)
    while len(newlist)>0: #since newlist is a list of lists we pop it until it is len==0
        y_item = (newlist.pop(0)) #pop the first list off
        y_item_tab = '\t'.join(y_item) #tab separate it
        fout.write(y_item_tab + '\n') #and write it to the file
fout.close() #close that file.