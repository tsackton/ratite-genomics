# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 12:24:50 2015

@author: phil-grayson
"""

import sys
import csv
phile = str(sys.argv[1])

class Vividict(dict): #used for autovivification of dictionary 
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value
        
cdic = Vividict() #make dictionary
fout = open('/Users/Phil/Desktop/maker/'+phile.split('.')[0]+'out.txt', 'w') #outfile
with open('/Users/Phil/Desktop/maker/'+phile, 'rU') as handle: #infile
    reader=csv.reader(handle,delimiter='\t') 
    species = str(handle).split('/')[-1].split('.')[0] #species name from infile 
    fout.write('GeneID'+'\t'+species+'_Lift'+'\t'+species+'_Inter'+'\n') #write the top line of outfile    
    for strLine in reader: 
        metaData = strLine[3] #metaData contains info from chicken - gene ID, genebank ID, start and stop of exon
        intersectCount = strLine[12] #intersectCount is the number of intersections between the pslbed and the RNAbed files
        GeneID = metaData.split(',')[0] #GeneID pops just the gene ID out of the metaData
        if GeneID not in cdic: #if the GeneID is not in cdic 
            cdic[GeneID][metaData]=intersectCount #add the ID, the metaData and the count
        else: #if the GeneID is in cdic
            if metaData not in cdic[GeneID]: #if the metaData is not in cdic
                cdic[GeneID][metaData]=intersectCount #add everything as above
            else: #if the metaData is in cdic
                if cdic[GeneID][metaData] == '0': #if the count is 0
                    cdic[GeneID][metaData]=intersectCount #add the current line's count - this will replace 0 with 0 or a higher value, but will not replace a higher value.  It provides us with a simple binary, 0 or not.
    for geneID in cdic: #once the whole file is read into the cdic, look at each geneID
        total=0 #start a total count
        positive=0 #and a count for intersectCount values above 0
        for exon in cdic[geneID]: #for each entry (exon) for a given geneID
            total+=1 #add 1 to total - this represents the number of exons lifted over from chicken
            if int(cdic[geneID][exon])>0: #and if the intersectCount is > 0
                positive+=1 #add 1 to positive as well
        fout.write(geneID.split(':')[1]+'\t'+str(total)+'\t'+str(positive)+'\n') #once all exons have been counted, write out the ID, the total, and the postive counts.
fout.close()
